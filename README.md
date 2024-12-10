# PINGVerter
A Python-based convertor for sonar logs collected with consumer-grade fishfinders.

`PINGVerter` is designed to pre-process sonar logs for [PINGMapper](https://cameronbodine.github.io/PINGMapper/) post-processing. It is not intended

## Compatibility

- [Humminbird](https://humminbird.johnsonoutdoors.com/) ( *.DAT/ *.SON/ *.IDX)
- [Lowrance](https://www.lowrance.com/) ( *.sl2 / *.sl3)

## Installation


```bash
pip install pingverter
```

## Usage - Example Only

`PINGVerter` *can* be used as a stand-alone engine for decoding sonar logs and exporting ping attributes to CSV. While this functionality is described here, please use [PINGMapper](https://cameronbodine.github.io/PINGMapper/) for all sonar log processing.

### Humminbird
```python
# Import
from pingverter import hum2pingmapper

# Parameters
inFile = r'C:\Path\To\Recording\Rec00001.DAT'
projDir = r'C:\Path\To\Outputs\MyProject'

sonar_object = hum2pingmapper(inFile, projDir)
```

### Lowrance
```python
# Import
from pingverter import low2pingmapper

# Parameters
inFile = r'C:\Path\To\Recording\Rec00001.DAT'
projDir = r'C:\Path\To\Outputs\MyProject'

sonar_object = low2pingmapper(inFile, projDir)
```

Outputs from the above examples will be exported to `C:\Path\To\Outputs\MyProject\meta`.

## Acknowledgments

`PINGVerter` has been made possible through mentorship, partnerships, financial support, open-source software, manuscripts, and documentation linked below.

*NOTE: The contents of this repository are those of the author(s) and do not necessarily represent the views of the individuals and organizations specifically mentioned here.*

- [Dr. Arthur Trembanis](https://www.udel.edu/academics/colleges/ceoe/departments/smsp/faculty/arthur-trembanis/) (Post-Doc Advisor) & [Dr. Daniel Buscombe](https://github.com/dbuscombe-usgs) (PhD Advisor)
- [Coastal Sediments, Hydrodynamics and Engineering Lab (CSHEL)](https://sites.udel.edu/ceoe-art/), [College of Earth, Ocean, & Environment (CEOE)](https://www.udel.edu/ceoe/), [University of Delaware](https://www.udel.edu/)
- [PINGMapper](https://cameronbodine.github.io/PINGMapper/) - Bodine, C. S., Buscombe, D., Best, R. J., Redner, J. A., & Kaeser, A. J. (2022). PING-Mapper: Open-source software for automated benthic imaging and mapping using recreation-grade sonar. Earth and Space Science, 9, e2022EA002469. https://doi.org/10.1029/2022EA002469
- [SL3Reader](https://github.com/halmaia/SL3Reader) - Halmai, Ákos; Gradwohl–Valkay, Alexandra; Czigány, Szabolcs; Ficsor, Johanna; Liptay, Zoltán Árpád; Kiss, Kinga; Lóczy, Dénes and Pirkhoffer, Ervin. 2020. "Applicability of a Recreational-Grade Interferometric Sonar for the Bathymetric Survey and Monitoring of the Drava River" ISPRS International Journal of Geo-Information 9, no. 3: 149. https://doi.org/10.3390/ijgi9030149
- [sonarlight](https://github.com/KennethTM/sonarlight) - Kenneth Thorø Martinsen
- [Navico (Lowrance, Simrad, B&G) Sonar Log File Format](https://www.memotech.franken.de/FileFormats/Navico_SLG_Format.pdf) - Herbert Oppmann
- [Vincent Capone](https://blacklaserlearning.com/) - Black Laser Learning


## Future Development, Collaborations, & Partnerships

If you are interested in partnering on future developments, please reach out to [Cameron Bodine, PhD](https://cameronbodine.github.io/).
