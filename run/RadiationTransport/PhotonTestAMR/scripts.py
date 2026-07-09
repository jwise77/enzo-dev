import yt
from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from yt.units.yt_array import YTQuantity, YTArray
# Calculate
# 1. Ionization front propagation and hydrodynamic response
# 2. Deviation from 1/r^2 in the photo-ionization field in the last output

first = 1
last = 25
NBINS=16
Myr = 3.1557e13
kpc = 3.086e21

########################################################################

def _MyRadius(field, data):
    center = data.get_field_parameter("center")
    dx = data["x"] - center[0]
    dy = data["y"] - center[1]
    dz = data["z"] - center[2]
    return np.sqrt(dx*dx + dy*dy + dz*dz)
yt.add_field(("gas", "Radius"), function=_MyRadius, take_log=False, units='', sampling_type="cell")

def _MyNeutralFrac(field, data):
    return data['H_p0_fraction'] / 0.75908798
yt.add_field(("gas", "Neutral_Fraction"), function=_MyNeutralFrac, take_log=False,
             units='', sampling_type="cell")
########################################################################

center = [1e-3]*3
time = []
radius = []

for i in range(first, last+1):
    amrfile = "DD%4.4d/data%4.4d" % (i,i)
    pf = yt.load(amrfile)

    x_bins_1d = 32
    r_min = pf.index.get_smallest_dx()
    r_max = pf.quan(1.0 - 1.0/64, 'code_length')
    
    sphere = pf.sphere(center, r_max)

    prof1d = yt.create_profile(sphere, ('index', 'radius'), fields=[("gas", "Neutral_Fraction"), ("enzo", "HI_kph")],
                               n_bins=x_bins_1d,
                               units = {('index', 'radius'):'code_length'})

    # Find the radius of the I-front (f_HI=0.5)
    res = np.abs(prof1d[("gas", "Neutral_Fraction")] - 0.5)
    ir = np.where(res == res.min())[0]
    r = np.interp(0.5, prof1d[("gas", "Neutral_Fraction")], prof1d.x.value)
    r = pf.quan(r, 'code_length')
    time.append(pf.current_time.to('s'))
    radius.append(r.to('cm'))

    del pf
    del prof1d

time = np.array(time)
radius = np.array(radius)

p = plt.subplot(111)
p.plot(time/Myr, radius/kpc, 'k-')
p.set_xlabel("Time (Myr)")
p.set_ylabel(r'$r_{\rm IF}$ (kpc)')
plt.savefig("IFrontRadius.png")

########################################################################
# Radial profiles for some outputs
########################################################################

all_profiles = {}
fields = [("gas", "density"), ("gas", "temperature"), ("gas", "H_p0_fraction"), ("gas", "H_p1_fraction")]
outputs = [3,5,10,15,25]
x_bins_1d = 20
r_min = 1.0/16
r_max = 1.0 - 1.0/16

for outp in outputs:
    amrfile = "DD%4.4d/data%4.4d" % (outp, outp)
    pf = yt.load(amrfile)
    sphere = pf.sphere(center, r_max)
    prof1d = yt.create_profile(sphere, ('index', 'radius'), fields, extrema={('index', 'radius'): (r_min,r_max)},
                               units = {('index', 'radius'):'code_length'}, n_bins=NBINS)                                   
    all_profiles[outp] = prof1d
    del pf

for f in fields:
    plt.clf()
    for outp in outputs:
        plt.semilogy(all_profiles[outp].x.value,
                     all_profiles[outp][f], label="%d Myr" % outp)
    plt.xlabel("Radius")
    plt.ylabel(f[1])
    plt.legend()
    plt.savefig("%sEvo.png" % f[1])
print("HII = ", all_profiles[25][("gas", "H_p1_fraction")])
########################################################################
# Some basic analysis on the final output
########################################################################

pf = yt.load("DD%4.4d/data%4.4d" % (last,last))
MyFields = [('enzo', 'HI_kph'), ('gas', 'H_p0_fraction'), ('gas', 'El_fraction'), ('gas', 'temperature'), ('gas', 'density')]
pc = yt.SlicePlot(pf,2,center=[0.5,0.5,1.0/64], fields=MyFields)

pc.save()
del pc

########################################################################
# Calculate deviation from 1/r^2 in the last output (inside 0.5*radius)
########################################################################

x_bins_1d = 20
r_min = 4*pf.index.get_smallest_dx()
r_max = 0.5*radius[-1]
r_max = YTQuantity(r_max, 'cm')
print("rmin, rmax = ", r_min, r_max)
sphere = pf.sphere(center, r_max)
fields = [('enzo', 'HI_kph')]
prof1d = yt.create_profile(sphere, ('index', 'radius'), fields, n_bins=NBINS,
                           extrema={('index', 'radius'): (r_min, r_max.to('code_length'))},
                           units = {('index', 'radius'): 'code_length'})
x_data = prof1d.x.value[2:-1]
y_data = prof1d[("enzo", "HI_kph")].value[2:-1]
valid = (x_data > 0) & (y_data > 0)
coeff, residual, tr1, tr2, tr3 = \
       np.polyfit(np.log(x_data[valid]),
                  np.log(y_data[valid]), 1, full=True)

print("="*72)
print("Inside 2*r_anyl: Radiation field slope = %f +/- %g" % \
      (coeff[0], float(residual[0]) if len(residual) > 0 else 0.0))
print("="*72)
