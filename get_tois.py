import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroquery.mast import Catalogs
from collections import OrderedDict

def is_daylight_savings(stop_date):
	year = int(stop_date[:4])
	month = int(stop_date[5:7])
	day = int(stop_date[8:])
	if (datetime(year, 3, 10) <= datetime(year, month, day) <= datetime(year, 10, 3)):
		return True
	else:
		return False

def gen_target_list(site, obs_type, start_date, start_time, stop_date, stop_time):
	# assign observing location and appropriate times
	if site == "lick":
		loc = EarthLocation(lat=37.3414*u.deg, lon=-121.6429*u.deg, height=1280*u.m)
		if is_daylight_savings(stop_date) is True:
			utcoffset = -7*60*u.minute
		else:
			utcoffset = -8*60*u.minute
	elif site == "maunakea":
		loc = EarthLocation(lat=19.8208*u.deg, lon=-155.4681*u.deg, height=4205*u.m)
		utcoffset = -10*60*u.minute

	midnight = Time(stop_date+" 00:00", format="iso", out_subfmt="date_hm") - utcoffset
	begin_time = Time(start_date+" "+start_time) - utcoffset
	end_time = Time(stop_date+" "+stop_time) - utcoffset
	delta_midnight = np.arange(int((begin_time-midnight).value*24*60), int((end_time-midnight).value*24*60), 1)*u.minute

	# find stars that will be observable
	TICIDs = []
	Kmags = []
	rmags = []
	start = []
	stop = []
	best = []
	am = []
	df = pd.read_csv("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe", sep="|")
	for i in tqdm(range(df.shape[0])):
		TOI = SkyCoord(df.RA[i], df.Dec[i], unit=(u.hourangle, u.deg))
		utcoffset = -7*60*u.minute 
		altazs = TOI.transform_to(AltAz(obstime=midnight+delta_midnight, location=loc)) 
		airmass = np.array(altazs.secz)
		best_times_idx = np.argwhere((airmass <= 2.5) & (airmass >= 1.0))[:,0]
		best_times = delta_midnight[best_times_idx]
		best_airmasses = airmass[best_times_idx]
		if best_times.shape[0] > 0:
			TICIDs.append(df["TIC ID"][i])
			star = Catalogs.query_object("TIC"+str(df["TIC ID"][i]), radius=1*u.arcsec, catalog="TIC").to_pandas()
			Kmags.append(star["Kmag"][0])
			rmags.append(star["rmag"][0])
			best_time_idx = np.argmin(best_airmasses)
			best_time = best_times[best_time_idx]
			best_airmass = best_airmasses[best_time_idx]  
			start.append((midnight+best_times[0]+utcoffset).value[11:])
			stop.append((midnight+best_times[-1]+utcoffset).value[11:])
			best.append((midnight+best_time+utcoffset).value[11:])
			am.append(np.round(best_airmass, 2))
	TICIDs = np.array(TICIDs)
	Kmags = np.array(Kmags)
	rmags = np.array(rmags)
	start = np.array(start)
	stop = np.array(stop)
	best = np.array(best)
	am = np.array(am)

	df2 = df[df["TIC ID"].isin(TICIDs)]
	df2 = df2.assign(**{"Kmag": Kmags, "rmag": rmags, "start time": start, "stop time": stop, "best time": best, "best airmass": am})
	df2.reset_index()
	new_RA = []
	new_Dec = []
	for i in range(len(np.array(df2['RA']))):
		new_RA.append(str(np.array(df2['RA'])[i]))
		new_Dec.append(str(np.array(df2['Dec'])[i]))

	# find previous observations of each observable star and combine everything into a new dataframe
	if obs_type == "spec":
		spec_obs = pd.read_csv("https://exofop.ipac.caltech.edu/tess/download_spect.php?sort=id&output=pipe", sep="|")
		observed = []
		telescope = []
		instrument = []
		user = []
		notes = []
		for i in range(len(TICIDs)):
			if TICIDs[i] in np.array(spec_obs["TIC ID"]):
				observed.append("Yes")
				idx = np.argwhere(TICIDs[i] == np.array(spec_obs["TIC ID"]))
				if len(idx) == 1:
					idx = idx[0,0]
					telescope.append(str(np.array(spec_obs["Telescope"])[idx]))
					instrument.append(str(np.array(spec_obs["Instrument"])[idx]))
					user.append(str(np.array(spec_obs["User"])[idx]))
					notes.append(str(np.array(spec_obs["Notes"])[idx]))
				elif len(idx) > 1:
					tel = ''
					inst = ''
					usr = ''
					nts = ''
					for j in range(len(idx)):
						tel += str(str(np.array(spec_obs["Telescope"])[idx[j]][0]))
						inst += str(str(np.array(spec_obs["Instrument"])[idx[j]][0]))
						usr += str(str(np.array(spec_obs["User"])[idx[j]][0]))
						nts += str(str(np.array(spec_obs["Notes"])[idx[j]][0]))  
						if j < len(idx)-1:
							tel += '|'
							inst += '|'
							itype += '|'
							usr += '|'
							nts += '|'  
					telescope.append(tel)
					instrument.append(inst)
					user.append(usr)
					notes.append(nts)           
			else:
				observed.append('')
				telescope.append('')
				instrument.append('')
				user.append('')
				notes.append('')

			d = OrderedDict(
				[("src", np.array(df2["Source"])),
				 ("TIC ID", np.array(df2["TIC ID"])),
				 ("TOI", np.array(df2["TOI"])),
				 ("TESS Disp", np.array(df2["TESS Disposition"])),
				 ("TFOPWG Disp", np.array(df2["TFOPWG Disposition"])),
				 ("SG3 Priority", np.array(df2["SG2"])),
				 ("RA", np.array(new_RA)),
				 ("Dec", np.array(new_Dec)),
				 ("T mag", np.array(df2["TESS Mag"])),
				 ("K mag", np.array(df2["Kmag"])),
				 ("r mag", np.array(df2["rmag"])),
				 ("start time", np.array(df2["start time"])),
				 ("stop time", np.array(df2["stop time"])),
				 ("best time", np.array(df2["best time"])),
				 ("best airmass", np.array(df2["best airmass"])),
				 ("Observed?", np.array(observed)),
				 ("Telescope", np.array(telescope)),
				 ("Instrument", np.array(instrument)),
				 ("Observer", np.array(user)),
				 ("Period (days)", np.array(df2["Period (days)"])),
				 ("Planet Radius (R_Earth)", np.array(df2["Planet Radius (R_Earth)"])),
				 ("Stellar Teff (K)", np.array(df2["Stellar Eff Temp (K)"])),
				 ("Stellar Radius (R_Sun)", np.array(df2["Stellar Radius (R_Sun)"])),
				 ("Sectors", np.array(df2["Sectors"])),
				 ("Created (UTC)", np.array(df2["Date TOI Created (UTC)"])),
				 ("Modified (UTC)", np.array(df2["Date TOI Modified (UTC)"])),
				 ("Comments", np.array(df2["Comments"]))  
				])

	elif obs_type == "img":
		img_obs = pd.read_csv("https://exofop.ipac.caltech.edu/tess/download_imaging.php?sort=id&output=pipe", sep="|")
		observed = []
		telescope = []
		instrument = []
		imgtype = []
		user = []
		notes = []
		for i in range(len(TICIDs)):
			if TICIDs[i] in np.array(img_obs["TIC ID"]):
				observed.append("Yes")
				idx = np.argwhere(TICIDs[i] == np.array(img_obs["TIC ID"]))
				if len(idx) == 1:
					telescope.append(str(np.array(img_obs["Telescope"])[idx[0,0]]))
					instrument.append(str(np.array(img_obs["Instrument"])[idx[0,0]]))
					imgtype.append(str(np.array(img_obs["Image Type"])[idx[0,0]]))
					user.append(str(np.array(img_obs["User"])[idx[0,0]]))
					notes.append(str(np.array(img_obs["Notes"])[idx[0,0]]))
				elif len(idx) > 1:
					tel = ''
					inst = ''
					itype = ''
					usr = ''
					nts = ''
					for j in range(len(idx)):
						tel += str(str(np.array(img_obs["Telescope"])[idx[j,0]]))
						inst += str(str(np.array(img_obs["Instrument"])[idx[j,0]]))
						itype += str(str(np.array(img_obs["Image Type"])[idx[j,0]]))
						usr += str(str(np.array(img_obs["User"])[idx[j,0]]))
						nts += str(str(np.array(img_obs["Notes"])[idx[j,0]]))
						if j < len(idx)-1:
							tel += '|'
							inst += '|'
							itype += '|'
							usr += '|'
							nts += '|' 
					telescope.append(tel)
					instrument.append(inst)
					imgtype.append(itype)
					user.append(usr)
					notes.append(nts)           
			else:
				observed.append('')
				telescope.append('')
				instrument.append('')
				imgtype.append('')
				user.append('')
				notes.append('')

		d = OrderedDict(
			[("src", np.array(df2["Source"])),
			 ("TIC ID", np.array(df2["TIC ID"])),
			 ("TOI", np.array(df2["TOI"])),
			 ("TESS Disp", np.array(df2["TESS Disposition"])),
			 ("TFOPWG Disp", np.array(df2["TFOPWG Disposition"])),
			 ("SG3 Priority", np.array(df2["SG3"])),
			 ("RA", np.array(new_RA)),
			 ("Dec", np.array(new_Dec)),
			 ("T mag", np.array(df2["TESS Mag"])),
			 ("K mag", np.array(df2["Kmag"])),
			 ("r mag", np.array(df2["rmag"])),
			 ("start time", np.array(df2["start time"])),
			 ("stop time", np.array(df2["stop time"])),
			 ("best time", np.array(df2["best time"])),
			 ("best airmass", np.array(df2["best airmass"])),
			 ("Observed?", np.array(observed)),
			 ("Telescope", np.array(telescope)),
			 ("Instrument", np.array(instrument)),
			 ("Image Type", np.array(imgtype)),
			 ("Observer", np.array(user)),
			 ("Period (days)", np.array(df2["Period (days)"])),
			 ("Planet Radius (R_Earth)", np.array(df2["Planet Radius (R_Earth)"])),
			 ("Stellar Teff (K)", np.array(df2["Stellar Eff Temp (K)"])),
			 ("Stellar Radius (R_Sun)", np.array(df2["Stellar Radius (R_Sun)"])),
			 ("Sectors", np.array(df2["Sectors"])),
			 ("Created (UTC)", np.array(df2["Date TOI Created (UTC)"])),
			 ("Modified (UTC)", np.array(df2["Date TOI Modified (UTC)"])),
			 ("Comments", np.array(df2["Comments"]))  
			])

	# save dataframe as a csv
	df3 = pd.DataFrame(data=d)
	df3.to_csv('TESS_targets.csv', index=False)

	return


if __name__ == "__main__":
	site = str(sys.argv[1])
	obs_type = str(sys.argv[2])
	start_date = str(sys.argv[3])
	start_time = str(sys.argv[4])
	stop_date = str(sys.argv[5])
	stop_time = str(sys.argv[6])
	gen_target_list(site, obs_type, start_date, start_time, stop_date, stop_time)