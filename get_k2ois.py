import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from collections import OrderedDict

def is_daylight_savings(stop_date):
	year = int(stop_date[:4])
	month = int(stop_date[5:7])
	day = int(stop_date[8:])
	if (datetime(year, 3, 8) <= datetime(year, month, day) <= datetime(year, 11, 1)):
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
	EPICIDs = []
	start = []
	stop = []
	best = []
	am = []
	df = pd.read_csv('https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All', skiprows=11)
	df = df[df['Poss_pc'] == 'Y']
	for i in tqdm(range(df.shape[0])):
		K2OI = SkyCoord(df.RA.values[i], df.Dec.values[i], unit=(u.hourangle, u.deg))
		utcoffset = -7*60*u.minute 
		altazs = K2OI.transform_to(AltAz(obstime=midnight+delta_midnight, location=loc)) 
		airmass = np.array(altazs.secz)
		best_times_idx = np.argwhere((airmass <= 2.5) & (airmass >= 1.0))[:,0]
		best_times = delta_midnight[best_times_idx]
		best_airmasses = airmass[best_times_idx]
		if best_times.shape[0] > 0:
			EPICIDs.append(df["EPIC ID"].values[i])
			best_time_idx = np.argmin(best_airmasses)
			best_time = best_times[best_time_idx]
			best_airmass = best_airmasses[best_time_idx]  
			start.append((midnight+best_times[0]+utcoffset).value[11:])
			stop.append((midnight+best_times[-1]+utcoffset).value[11:])
			best.append((midnight+best_time+utcoffset).value[11:])
			am.append(np.round(best_airmass, 2))
	EPICIDs = np.array(EPICIDs)
	start = np.array(start)
	stop = np.array(stop)
	best = np.array(best)
	am = np.array(am)

	df2 = df[df["EPIC ID"].isin(EPICIDs)]
	df2 = df2.assign(**{"start time": start, "stop time": stop, "best time": best, "best airmass": am})
	df2.reset_index()
	new_RA = []
	new_Dec = []
	for i in range(len(np.array(df2['RA']))):
		new_RA.append(str(np.array(df2['RA'])[i]))
		new_Dec.append(str(np.array(df2['Dec'])[i]))

	# find previous observations of each observable star and combine everything into a new dataframe
	if obs_type == "spec":
		spec_obs = pd.read_csv("https://exofop.ipac.caltech.edu/k2/download_spect.php?sort=id&output=csv", usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13])
		observed = []
		telescope = []
		instrument = []
		user = []
		notes = []
		for i in range(len(EPICIDs)):
			if EPICIDs[i] in np.array(spec_obs["EPIC ID"]):
				observed.append("Yes")
				idx = np.argwhere(EPICIDs[i] == np.array(spec_obs["EPIC ID"]))
				if len(idx) == 1:
					telescope.append(str(np.array(spec_obs["Telescope"])[idx[0,0]]))
					instrument.append(str(np.array(spec_obs["Instrument"])[idx[0,0]]))
					user.append(str(np.array(spec_obs["User"])[idx[0,0]]))
					notes.append(str(np.array(spec_obs["Notes"])[idx[0,0]]))
				elif len(idx) > 1:
					tel = ''
					inst = ''
					itype = ''
					usr = ''
					nts = ''
					for j in range(len(idx)):
						tel += str(str(np.array(spec_obs["Telescope"])[idx[j,0]]))
						inst += str(str(np.array(spec_obs["Instrument"])[idx[j,0]]))
						usr += str(str(np.array(spec_obs["User"])[idx[j,0]]))
						nts += str(str(np.array(spec_obs["Notes"])[idx[j,0]]))    
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
			[("EPIC ID", np.array(df2["EPIC ID"])),
			 ("Num_CP", np.array(df2["Num_plan"])),
			 ("RA", np.array(new_RA)),
			 ("Dec", np.array(new_Dec)),
			 ("Kep mag", np.array(df2["Mag_Kep"])),
			 ("K mag", np.array(df2["Mag_Ks"])),
			 ("start time", np.array(df2["start time"])),
			 ("stop time", np.array(df2["stop time"])),
			 ("best time", np.array(df2["best time"])),
			 ("best airmass", np.array(df2["best airmass"])),
			 ("lowest zd", np.arccos(1/np.array(df2["best airmass"])) * 180/np.pi),
			 ("Observed?", np.array(observed)),
			 ("Telescope", np.array(telescope)),
			 ("Instrument", np.array(instrument)),
			 ("Observer", np.array(user)),
			 ("Last Modified", np.array(df2["Last Mod"])), 
			])

	elif obs_type == "img":
		img_obs = pd.read_csv("https://exofop.ipac.caltech.edu/k2/download_imaging.php?sort=id&output=csv", usecols=[0,1,2,3,4,5,6,7,8,9,10])
		observed = []
		telescope = []
		instrument = []
		user = []
		notes = []
		for i in range(len(EPICIDs)):
			if EPICIDs[i] in np.array(img_obs["EPIC ID"]):
				observed.append("Yes")
				idx = np.argwhere(EPICIDs[i] == np.array(img_obs["EPIC ID"]))
				if len(idx) == 1:
					telescope.append(str(np.array(img_obs["Telescope"])[idx[0,0]]))
					instrument.append(str(np.array(img_obs["Instrument"])[idx[0,0]]))
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
					user.append(usr)
					notes.append(nts)  
			else:
				observed.append('')
				telescope.append('')
				instrument.append('')
				user.append('')
				notes.append('')

		d = OrderedDict(
			[("EPIC ID", np.array(df2["EPIC ID"])),
			 ("Num_CP", np.array(df2["Num_plan"])),
			 ("RA", np.array(new_RA)),
			 ("Dec", np.array(new_Dec)),
			 ("Kep mag", np.array(df2["Mag_Kep"])),
			 ("K mag", np.array(df2["Mag_Ks"])),
			 ("start time", np.array(df2["start time"])),
			 ("stop time", np.array(df2["stop time"])),
			 ("best time", np.array(df2["best time"])),
			 ("best airmass", np.array(df2["best airmass"])),
			 ("lowest zd", np.arccos(1/np.array(df2["best airmass"])) * 180/np.pi),
			 ("Observed?", np.array(observed)),
			 ("Telescope", np.array(telescope)),
			 ("Instrument", np.array(instrument)),
			 ("Observer", np.array(user)),
			 ("Last Modified", np.array(df2["Last Mod"])), 
			])

	# save dataframe as a csv
	df3 = pd.DataFrame(data=d)
	df3.to_csv('K2_targets.csv', index=False)

	return


if __name__ == "__main__":
	site = str(sys.argv[1])
	obs_type = str(sys.argv[2])
	start_date = str(sys.argv[3])
	start_time = str(sys.argv[4])
	stop_date = str(sys.argv[5])
	stop_time = str(sys.argv[6])
	gen_target_list(site, obs_type, start_date, start_time, stop_date, stop_time)