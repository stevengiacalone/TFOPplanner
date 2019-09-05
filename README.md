# TFOPplanner
Tools that assist with the selection of observable TESS Objects of Interest for ground based follow-up. 
___________________________________
### get_tois.py / get_k2ois.py
Pulls data from exofop to generate a CSV of objects of interest that will be observable from a specified observatory within a specified time interval. CSV includes useful planning info, such as SG priority and previous observations of each object.

Run with: 
'''
python3 get_tois.py SITE OBS_TYPE START_DATE START_TIME END_DATE END_TIME
'''
where SITE is the observatory (lick or maunakea), OBS_TYPE is the type of observation being taken (spec or img), START_DATE/END_DATE is the date that you will start/end you observations in the timezone of the observatory (formatted as YYYY-MM-DD), and START_TIME/END_TIME is the time that you will start/end your observations in the timezone of the observatory (formatted as HH:MM).

### to_jskycalc.py
Converts a CSV with a list of object names and coordinates into a format readable by jskycalc.

Run with:
'''
python3 to_jskycalc.py INPUT.CSV
'''
where INPUT.CSV is a CSV with the same headers as those in files created with the above scripts.
