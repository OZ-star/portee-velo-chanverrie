# À portée de vélo – Chanverrie

Carte interactive des temps de trajet à vélo sur la commune de Chanverrie (Vendée, code INSEE 85302).

Déplacez la souris sur la carte : les temps de trajet depuis le curseur se recalculent en direct. Cliquez pour fixer le départ, puis cliquez ailleurs pour obtenir l'itinéraire et sa durée. Sur mobile, faites glisser le repère D.

Adaptation du code d'[À portée de vélo](https://github.com/sit-pmo/portee-velo) (© PETR Pôle Marennes-Oléron, licence MIT), lui-même inspiré d'[À portée de tram](https://tram.camilleroux.com/) de Camille Roux.

## Ce que fait la carte

- Mode **Pistes + routes** : le réseau aménagé complété par les routes et chemins, avec un temps majoré selon le type de voie.
- Mode **Pistes et bandes** : uniquement le réseau cyclable aménagé. Les zones grises sont à plus de 800 m de ce réseau.
- Le calcul déborde sur les communes voisines (jusqu'à 3 km, mais seulement 0,6 km au sud et 1,5 km au nord, limites de l'extrait OpenStreetMap utilisé), pour ne pas couper les itinéraires en limite. Les pourcentages de surface atteinte portent sur la commune seule.

## Modèle de calcul

| Élément                              | Hypothèse                       |
| ------------------------------------ | ------------------------------- |
| Piste cyclable, route, rue           | 20 km/h                         |
| Chemin (`track`)                     | 15 km/h, temps majoré de 50 %   |
| Sentier (`path`)                     | 13 km/h, temps majoré de 50 %   |
| Route ou rue partagée                | temps majoré de 50 %            |
| Route secondaire                     | temps doublé                    |
| Route principale sans aménagement    | temps triplé                    |
| Sens unique pris à contre-sens       | vélo poussé, à 5 km/h           |
| Accès au réseau depuis un point      | à vol d'oiseau, majoré de 30 %, à 16 km/h |

Ni vent ni relief. Exclus : autoroutes et bretelles, chantiers, escaliers, cheminements piétons.

## Limites des données

Les données viennent uniquement d'OpenStreetMap (extrait BBBike au format shapefile). Cet extrait ne contient ni le revêtement ni les bandes cyclables tracées sur les routes : seules les voies de type `cycleway` sont comptées comme pistes. Le mode « Pistes et bandes » est donc très partiel tant qu'une base cyclable locale n'est pas intégrée.

## Fabrication des données

`build_graph_velo.py` lit `roads.shp` (BBBike) et le contour communal (`contour.geojson`, geo.api.gouv.fr) et produit `data/portee_velo.js` : graphe du réseau et grille de cellules de 100 m. Dépendances : `numpy`, `scipy`.

```
python3 build_graph_velo.py --roads osm/roads --contour contour.geojson --out data/portee_velo.js
```

La page est statique (`index.html`, `coeur.js`, `data/portee_velo.js`) : elle peut être servie par GitHub Pages sans serveur ni base de données.

## Sources et licences

- Routes, chemins et pistes cyclables : © contributeurs [OpenStreetMap](https://www.openstreetmap.org/copyright) (ODbL). Les données dérivées contenues dans `data/portee_velo.js` sont diffusées sous ODbL.
- Contour communal : [geo.api.gouv.fr](https://geo.api.gouv.fr).
- Fond de carte et recherche d'adresses : [Géoplateforme IGN](https://geoservices.ign.fr).
- Bibliothèques : [Leaflet](https://leafletjs.com), [d3](https://d3js.org).

Le code (`index.html`, `coeur.js`, `build_graph_velo.py`) est diffusé sous [licence MIT](LICENSE). Les données ne sont pas couvertes par cette licence.
