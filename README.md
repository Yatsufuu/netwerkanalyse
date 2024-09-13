# Rijnland - Netwerk analyse

## Script werkwijze

Dit script analyseert het verloop van de gemengde stelsel / vuilwaterriool om te bepalen welke percelen waarschijnlijk aangesloten zijn op de riolering. 

De werkwijze is als volgt;

* ophalen data uit online services (met name BAG)
* automatisch bepalen welke leidingen bij welke pompen horen
* verwijderen locaties waar de leiding een waterpartij oversteekt
* creeeren buffer rondom op een pomp aangesloten leidingen
* bepaling raakvlak tussen buffers en percelen

Het ophalen van de data gebeurt via het ```download_data.py``` script. De gemeentes waarvoor de data gedownload wordt kunnen gedefinieerd worden via het ```settings.py``` bestand. Het script maakt zelf de mappen aan om de data in weg te schrijven.

**LET OP** dit script maakt gebruik van online services die met name onder de BAG vallen. Deze zijn aan veranderingen onderhevig waardoor het op zijn tijd nodig is om te kijken of de services aangepast zijn. Dit vereist vaak het wijzigen van URL's en soms veldnamen en kan tijdrovend zijn.

De analyse vindt plaats via het ```analyse.py``` script. In dit script staat in de code aangegeven welke stappen worden doorlopen. Voor iedere stap wordt in de uitvoer map een GIS bestand gemaakt dat in QGis of ArcGIS geopened kan worden om te kijken wat voor effect de stap heeft gehad. Op deze manier kan gekeken worden of bepaalde stappen verder geoptimaliseerd kunnen worden en kan eventueel debugging plaatsvinden als er rare gegevens worden gegenereerd. 

De uiteindelijke uitvoerlaag wordt in de analyse map gegenereerd en heet ```05_plots_with_sewer```. In deze laag zijn de volgens de analyse aangesloten percelen weergegeven.

## Instellingen

In overleg met Rijnland zijn de volgende belangrijke uitgangspunten aangehouden;

* leidingen kunnen geen waterpartijen oversteken
* de afstand tussen een pomp en een leiding bedraagt maximaal 10 meter (1)
* de afstand tussen twee aaneengesloten leidingen bedraagt maximaal 1 meter (2)
* de afstand tussen een aangesloten perceel en de leiding bedraagt maximaal 40 meter (3)

ad 1) een pomp heeft een x,y coordinaat, de leiding heeft een x,y coordinaat voor het begin en het eind, er wordt gezocht naar de kortste afstand tussen pomp en begin of einde leiding
ad 2) in de data zijn de begin- en eindpunten van twee leidingen niet altijd precies gelijk waardoor het nodig is om een maximale afstand aan te houden
ad 3) de buffer rondom de leidingen wordt als aangesloten op een pomp beschouwd als er in de buffer een actieve (op leidingen aangesloten) pomp aanwezig is

De instellingen zijn te vinden in het bestand ```settings.py``` en kunnen aangepast worden. 

**LET OP** Indien de instellingen gewijzigd zijn is het nodig om in het ```analyse.py``` script de waarde voor ```FORCE_RELOAD``` op True te zetten. Deze waarde zorgt er in eerste instantie voor dat het zware werk zoveel mogelijk gecached wordt zodat het rekentijd scheelt. Bij nieuwe instellingen moet de cache echter overschreven worden omdat de uitgangspunten gewijzigd zijn. Als de analyse eenmaal klaar is kan ```FORCE_RELOAD``` weer op False worden gezet zolang de uitgangspunten niet wijzigen.

## Installatie

Bij het gebruik van git gebruik;

```git clone git@github.com:breinbaas/netwerkanalyse.git``` 

Of download het ZIP bestand via;

https://github.com/breinbaas/netwerkanalyse | Code | Download ZIP

Open de map en open een console in de map

Creeer (in de werkmap) een virtuele Python omgeving;

```python -m venv .venv```

Activeer de omgeving

```.venv/Scripts/activate```

En installeer de vereiste packages

```python -m pip install -r requirements.txt```

**LET OP** de GDAL package zorgt vrijwel altijd voor problemen op Windows omgevingen. Dit komt omdat er code gecompileerd moet worden en daar zijn afhankelijkheden voor nodig die lang niet altijd op de Windows omgeving staan. Als oplossing daarvoor kan beter het voorgecompileerde bestand worden gebruikt. Ga hiervoor naar;

https://github.com/cgohlke/geospatial-wheels/releases

en download de versie van GDAL die bij je Python versie hoort. 

Je kunt je Python versie bepalen door ```python``` in de commandline uit te voeren. In de regels die verschijnen zie je de huidige versie (bv 3.12.4). Gebruik ```exit()``` om weer uit de Python shell te komen.

Zoek binnen de eerder genoemde link naar GDAL-3.8.4-cp{PYTHON VERSIE} en dan de win_amd64.whl bijvoorbeeld;

GDAL-3.8.4-cp312-cp312-win_amd64.whl

Voor Python versie 3.12

Let op dat we er hier gemakshalve vanuitgaan dat er op een 64bit Windows versie wordt gewerkt wat al lange tijd standaard is. Mocht er een 32 bits systeem worden gebruikt (**niet aannemelijk**) dan is de win32 versie nodig.

Installleer dit pakket eenvoudig via;

```python -m pip install {PAD NAAR WHL BESTAND}```

bv

```python -m pip install C:\\User\\ItsMe\\Downloads\\GDAL-3.8.4-cp312-cp312-win_amd64.whl```

Zolang de scripts uitgevoerd worden binnen de aangemaakte virtuele omgeving moeten ze over alle benodigde packages kunnen beschikken.

## Opmerkingen

Dit script is `as is` en expert judgement is nodig om het script en de uitkomsten te kunnen interpreteren.






