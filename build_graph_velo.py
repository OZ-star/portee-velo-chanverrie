#!/usr/bin/env python3
"""Prépare data/portee_velo.js pour « À portée de vélo – Chanverrie ».

Entrées  : extrait OpenStreetMap BBBike (shapefile roads.*) + contour communal (GeoJSON, geo.api.gouv.fr).
Sortie   : data/portee_velo.js (même format que le dépôt sit-pmo/portee-velo, lu par coeur.js).

Usage    : python3 build_graph_velo.py --roads osm/roads --contour contour.geojson --out out/data/portee_velo.js
Sans --contour, un cercle d'essai est utilisé : sortie réservée aux tests, jamais à publier.
Dépendances : numpy, scipy (aucune bibliothèque géospatiale).
"""
import argparse, collections, datetime, json, math, os, sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shp

R_TERRE = 6371008.8
RAD = math.pi / 180
M_PAR_DEG = 111319.49          # même valeur que le dépôt d'origine (pas de 100 m)
PAS_M = 100
MARGE_M = 3000                 # le calcul déborde de 3 km autour du territoire
RACCORD_M = 15                 # extrémités de piste reliées à la voie la plus proche sous 15 m

# Classes : 0 piste, 1 bande, 2 route/chemin, 3 secondaire, 4 principale
# type OSM -> (classe, vitesse km/h). Tout ce qui n'est pas listé est exclu
# (autoroutes, bretelles, chantiers, escaliers, cheminements piétons).
TYPES = {
    "cycleway":       (0, 20),
    "residential":    (2, 20), "unclassified": (2, 20), "living_street": (2, 20),
    "service":        (2, 20), "tertiary": (2, 20), "tertiary_link": (2, 20), "road": (2, 20),
    "track":          (2, 15), "path": (2, 13),
    "secondary":      (3, 20), "secondary_link": (3, 20),
    "primary":        (4, 20), "primary_link": (4, 20), "trunk": (4, 20), "trunk_link": (4, 20),
}
ORIGINES = ["osm", "raccord"]


def dist_m(lon1, lat1, lon2, lat2):
    """Même formule que distanceM() de coeur.js."""
    x = (lon2 - lon1) * RAD * np.cos((lat1 + lat2) / 2 * RAD)
    y = (lat2 - lat1) * RAD
    return np.hypot(x, y) * R_TERRE


def charger_contour(chemin, lon0, lat0):
    if chemin:
        with open(chemin, encoding="utf-8") as f:
            gj = json.load(f)
        geom = gj["geometry"] if gj.get("type") == "Feature" else (
            gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj)
    else:
        print("ATTENTION : pas de contour fourni, cercle d'essai de 4,3 km (NE PAS PUBLIER).", file=sys.stderr)
        pts = [[lon0 + 4300 / (M_PAR_DEG * math.cos(lat0 * RAD)) * math.cos(t),
                lat0 + 4300 / M_PAR_DEG * math.sin(t)] for t in np.linspace(0, 2 * math.pi, 181)]
        geom = {"type": "Polygon", "coordinates": [pts]}
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    return {"type": "MultiPolygon", "coordinates": polys}


def anneaux(contour):
    return [np.array(r, dtype=float) for poly in contour["coordinates"] for r in poly]


def dans_polygone(lon, lat, contour):
    """Remplissage pair-impair (trous gérés) sur l'ensemble des anneaux d'un MultiPolygon."""
    dedans = np.zeros(lon.shape, dtype=bool)
    for poly in contour["coordinates"]:
        for r in poly:
            r = np.array(r, dtype=float)
            x1, y1, x2, y2 = r[:-1, 0], r[:-1, 1], r[1:, 0], r[1:, 1]
            for i in range(len(x1)):
                cond = (y1[i] > lat) != (y2[i] > lat)
                if not cond.any():
                    continue
                xi = x1[i] + (lat - y1[i]) * (x2[i] - x1[i]) / (y2[i] - y1[i] + 1e-300)
                dedans ^= cond & (lon < xi)
    return dedans


def dist_au_bord(lon, lat, contour, lon0, lat0):
    """Distance (m) de chaque point au contour (segments), par paquets."""
    kx = RAD * R_TERRE * math.cos(lat0 * RAD)
    ky = RAD * R_TERRE
    px, py = (lon - lon0) * kx, (lat - lat0) * ky
    best = np.full(lon.shape, np.inf)
    for r in anneaux(contour):
        ax, ay = (r[:-1, 0] - lon0) * kx, (r[:-1, 1] - lat0) * ky
        bx, by = (r[1:, 0] - lon0) * kx, (r[1:, 1] - lat0) * ky
        dx, dy = bx - ax, by - ay
        l2 = dx * dx + dy * dy + 1e-12
        for s in range(0, len(ax), 40):
            sl = slice(s, s + 40)
            t = ((px[:, None] - ax[sl]) * dx[sl] + (py[:, None] - ay[sl]) * dy[sl]) / l2[sl]
            t = np.clip(t, 0, 1)
            qx, qy = ax[sl] + t * dx[sl], ay[sl] + t * dy[sl]
            d = np.hypot(px[:, None] - qx, py[:, None] - qy).min(axis=1)
            best = np.minimum(best, d)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roads", default="osm/roads")
    ap.add_argument("--contour")
    ap.add_argument("--out", default="out/data/portee_velo.js")
    ap.add_argument("--centre", default="-0.9857,46.9468", help="lon,lat du départ par défaut")
    ap.add_argument("--date-osm", default=datetime.date.today().isoformat())
    a = ap.parse_args()

    lon0, lat0 = map(float, a.centre.split(","))
    geoms = shp.read_shp(a.roads + ".shp")
    attrs = shp.read_dbf(a.roads + ".dbf")

    # ---------- 1. Arêtes brutes, nœuds uniques (coordonnées exactes à 1e-7°) ----------
    cle2id, nlon, nlat = {}, [], []
    def noeud(x, y):
        k = (round(x * 1e7), round(y * 1e7))
        i = cle2id.get(k)
        if i is None:
            i = cle2id[k] = len(nlon); nlon.append(x); nlat.append(y)
        return i

    noms, nom2id = [""], {"": 0}
    def id_nom(n):
        if n not in nom2id:
            nom2id[n] = len(noms); noms.append(n)
        return nom2id[n]

    A, B, CLS, SENS, VIT, NOM, ORI, EXTR = [], [], [], [], [], [], [], []
    exclus = collections.Counter()
    for g, r in zip(geoms, attrs):
        if g is None or r is None:
            continue
        t = r["type"]
        if t not in TYPES:
            exclus[t] += 1
            continue
        cls, vit = TYPES[t]
        sens = 1 if r.get("oneway") == 1 else (2 if r.get("oneway") == -1 else 0)
        nom = id_nom(r.get("name") or r.get("ref") or "")
        for part in g["parts"]:
            ids = [noeud(x, y) for x, y in part]
            for u, v in zip(ids[:-1], ids[1:]):
                if u == v:
                    continue
                A.append(u); B.append(v); CLS.append(cls); SENS.append(sens); VIT.append(vit)
                NOM.append(nom); ORI.append(0)
            if cls == 0 and len(ids) > 1:
                EXTR.append((ids[0], ids[-1], r["osm_id"]))
    nlon, nlat = np.array(nlon), np.array(nlat)
    print("arêtes brutes :", len(A), "| nœuds :", len(nlon), "| types exclus :", dict(exclus))

    # ---------- 2. Raccords : extrémités de piste -> nœud de voie le plus proche (< 15 m) ----------
    kx = RAD * R_TERRE * math.cos(lat0 * RAD); ky = RAD * R_TERRE
    xy = np.column_stack([(nlon - lon0) * kx, (nlat - lat0) * ky])
    deg = np.bincount(np.array(A + B), minlength=len(nlon))
    est_route = np.zeros(len(nlon), bool)
    for u, v, c in zip(A, B, CLS):
        if c >= 2:
            est_route[u] = est_route[v] = True
    arbre_routes = cKDTree(xy[est_route]); ids_routes = np.where(est_route)[0]
    voisins = collections.defaultdict(set)
    for u, v in zip(A, B):
        voisins[u].add(v); voisins[v].add(u)
    nr = 0
    for u in {x for e in EXTR for x in e[:2]}:
        if deg[u] != 1 or est_route[u]:
            continue
        d, k = arbre_routes.query(xy[u], k=3)
        for dd, kk in zip(np.atleast_1d(d), np.atleast_1d(k)):
            if dd <= RACCORD_M and kk < len(ids_routes) and ids_routes[kk] != u and ids_routes[kk] not in voisins[u]:
                A.append(u); B.append(int(ids_routes[kk])); CLS.append(0); SENS.append(0); VIT.append(20)
                NOM.append(0); ORI.append(1); nr += 1
                break
    print("raccords ajoutés :", nr)

    A, B = np.array(A), np.array(B)
    CLS, SENS, VIT, NOM, ORI = map(np.array, (CLS, SENS, VIT, NOM, ORI))

    # ---------- 3. Composante principale (les fragments isolés fausseraient l'accroche du départ) ----------
    n = len(nlon)
    mat = coo_matrix((np.ones(len(A)), (A, B)), shape=(n, n))
    ncomp, lab = connected_components(mat, directed=False)
    taille_arete = np.bincount(lab[A], minlength=ncomp)
    principale = int(np.argmax(taille_arete))
    print("composantes :", ncomp, "| arêtes de la principale :", int(taille_arete[principale]),
          "| autres (arêtes) :", sorted(taille_arete.tolist(), reverse=True)[1:8])
    garde = lab[A] == principale
    A, B, CLS, SENS, VIT, NOM, ORI = A[garde], B[garde], CLS[garde], SENS[garde], VIT[garde], NOM[garde], ORI[garde]
    used = np.unique(np.concatenate([A, B]))
    remap = -np.ones(n, dtype=np.int64); remap[used] = np.arange(len(used))
    nlon, nlat, xy = nlon[used], nlat[used], xy[used]
    A, B = remap[A], remap[B]
    # nœuds en 1e-5° comme le dépôt d'origine
    lon5 = np.round(nlon * 1e5).astype(int); lat5 = np.round(nlat * 1e5).astype(int)
    L_dm = np.round(dist_m(nlon[A], nlat[A], nlon[B], nlat[B]) * 10).astype(int)
    L_dm = np.maximum(L_dm, 1)

    km = [round(float(L_dm[CLS == c].sum()) / 10000, 1) for c in range(5)]
    print("km par classe (piste, bande, route, secondaire, principale) :", km)
    print("arêtes :", len(A), "| nœuds :", len(nlon))

    # ---------- 4. Contour et grille de 100 m ----------
    contour = charger_contour(a.contour, lon0, lat0)
    pts = np.vstack(anneaux(contour))
    marge_lon = MARGE_M / (M_PAR_DEG * math.cos(lat0 * RAD)); marge_lat = MARGE_M / M_PAR_DEG
    dlat = PAS_M / M_PAR_DEG; dlon = dlat / math.cos(lat0 * RAD)
    w = math.floor((pts[:, 0].min() - marge_lon - dlon) / dlon) * dlon
    top = math.ceil((pts[:, 1].max() + marge_lat + dlat) / dlat) * dlat
    east = pts[:, 0].max() + marge_lon + dlon; south = pts[:, 1].min() - marge_lat - dlat
    nx = int(math.ceil((east - w) / dlon)); ny = int(math.ceil((top - south) / dlat))
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny))
    clon = (w + (ii.ravel() + 0.5) * dlon); clat = (top - (jj.ravel() + 0.5) * dlat)
    cxy = np.column_stack([(clon - lon0) * kx, (clat - lat0) * ky])
    ile = dans_polygone(clon, clat, contour)
    bord = dist_au_bord(clon, clat, contour, lon0, lat0)
    dedans = ile | (bord <= MARGE_M)
    # La zone de calcul ne dépasse pas l'emprise de l'extrait routier : au-delà il n'y a pas de données,
    # et les cellules s'afficheraient « hors d'atteinte » à tort.
    bx = [p[0] for gg in geoms if gg for part in gg["parts"] for p in part]
    by = [p[1] for gg in geoms if gg for part in gg["parts"] for p in part]
    emprise = (clon >= min(bx)) & (clon <= max(bx)) & (clat >= min(by)) & (clat <= max(by))
    dedans &= emprise
    print("grille :", nx, "x", ny, "| cellules dedans :", int(dedans.sum()), "| sur le territoire :", int(ile.sum()),
          "| surface ≈ %.1f km²" % (ile.sum() * PAS_M * PAS_M / 1e6))

    def plus_proche(masque_noeuds):
        ids = np.where(masque_noeuds)[0]
        d, k = cKDTree(xy[ids]).query(cxy)
        return ids[k].astype(int), np.round(d).astype(int)

    actif_tout = np.zeros(len(nlon), bool); actif_tout[A] = True; actif_tout[B] = True
    actif_amen = np.zeros(len(nlon), bool)
    am = CLS <= 1; actif_amen[A[am]] = True; actif_amen[B[am]] = True
    pro_tout, dist_tout = plus_proche(actif_tout)
    if actif_amen.any():
        pro_amen, dist_amen = plus_proche(actif_amen)
    else:
        pro_amen = -np.ones(nx * ny, int); dist_amen = np.full(nx * ny, 99999)

    # ---------- 5. Écriture ----------
    contour_arr = {"type": "MultiPolygon",
                   "coordinates": [[[[round(x, 5), round(y, 5)] for x, y in r] for r in poly] for poly in contour["coordinates"]]}
    attr = (CLS + (SENS << 3) + (VIT << 5)).astype(int)
    D = {
        "meta": {"genere": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "osm": a.date_osm, "pas_m": PAS_M,
                 "depart_defaut": [lon0, lat0], "km_par_classe": km, "origines": ORIGINES,
                 "reseau": "OpenStreetMap (extrait BBBike)"},
        "contour": contour_arr,
        "noeuds": {"lon": lon5.tolist(), "lat": lat5.tolist()},
        "aretes": {"a": A.tolist(), "b": B.tolist(), "l": L_dm.tolist(), "attr": attr.tolist(),
                   "nom": NOM.tolist(), "o": ORI.tolist()},
        "noms": noms,
        "grille": {"w": w, "n": top, "dlon": dlon, "dlat": dlat, "nx": nx, "ny": ny,
                   "dedans": dedans.astype(int).tolist(), "ile": ile.astype(int).tolist(),
                   "pro_tout": pro_tout.tolist(), "dist_tout": dist_tout.tolist(),
                   "pro_amen": pro_amen.tolist(), "dist_amen": dist_amen.tolist()},
    }
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write("window.PORTEE_VELO=" + json.dumps(D, separators=(",", ":"), ensure_ascii=False) + ";")
    print("écrit :", a.out, "(%.1f Mo)" % (os.path.getsize(a.out) / 1e6))


if __name__ == "__main__":
    main()
