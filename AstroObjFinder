# Python Program: ASTRONOMY OBJECT FINDER #
# Author: Todd Chevrier - Amateur Astronomer - Citizen Scientist
# This program runs on the slow side when trying to interact with the graph
# Otherwise leave it alone and it will run for as long as it is open 
# and updates the graph every 5 seconds.
##################################################################
# Python standard-library
from urllib.parse import urlencode
from urllib.request import urlretrieve
import urllib.request
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Third-party dependencies **You may need to install them ** Install matplotlib,numpy, astropy, IPython, PIL **
from matplotlib.animation import FuncAnimation
from datetime import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from IPython.display import display, Image
from astropy.time import Time
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from PIL import Image
from astropy.visualization import astropy_mpl_style, quantity_support
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use(astropy_mpl_style)
quantity_support()

#Set my location and time zone **SET YOUR LOCATION** Change the x's
My_Location =  EarthLocation(lat=xx.xx*u.deg, lon=-xx.xx*u.deg, height=xxx*u.m)

#Set Time
utcoffset = -5*u.hour  # Central Daylight Time  **Change this if you are not in the US Central Time Zone
ctime = Time(datetime.now()) - utcoffset

# initialize a SkyCood object named obj_center
obj = input("Enter Object Name: ")
obj_center = SkyCoord.from_name(obj)
my_Loc_altaz = obj_center.transform_to(AltAz(obstime=ctime,location=My_Location))

#Get Object Image
url = "http://archive.stsci.edu/cgi-bin/dss_search?f=GIF&ra=" + obj_center.ra.__str__() + "&dec=" + obj_center.dec.__str__()
image = Image.open(urllib.request.urlopen(url))
width, height = image.size
image.show()

#Print the initial request    
print(obj + ' Location is: ')
print(obj_center.ra, obj_center.dec)
print(obj_center.ra.hour, obj_center.dec)
print('=====================')
print(" Altitude = {0.alt:.4}".format(my_Loc_altaz) + "\n Azimuth = {0.az:.4}".format(my_Loc_altaz))
print('=====================')
print("Image Location: " + url)
print(obj_center)

#Create the pretty graphs
fig, ax = plt.subplots()
fig.canvas.set_window_title('Astronomy Object Finder')

def animate(i):
    ax.clear()

    #Set Current Time
    ctime = Time(datetime.now()) - utcoffset
    obj_center = SkyCoord.from_name(obj)
    my_Loc_altaz = obj_center.transform_to(AltAz(obstime=ctime,location=My_Location))

    curpt = "    Alt = {0.alt:.5}".format(my_Loc_altaz) + "\n    Az = {0.az:.6}".format(my_Loc_altaz) + "\n <--"

    curAlt = ("{.alt:.6}".format(my_Loc_altaz))[:5]
    curAz =  ("{.az:.7}".format(my_Loc_altaz))[:6]

    midnight = Time(datetime.now()) - utcoffset
    delta_midnight = np.linspace(-2, 10, 100)*u.hour

    frame_Current = AltAz(obstime=midnight+delta_midnight, location=My_Location)

    m33altazs_Current = my_Loc_altaz.transform_to(frame_Current)

    m33airmasss_Current = m33altazs_Current.secz

    #Get Sun Position
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_Current = midnight + delta_midnight
    frame_Current = AltAz(obstime=times_Current, location=My_Location)
    sunaltazs_Current = get_sun(times_Current).transform_to(frame_Current)

    #Get Moon Position
    moon_Current = get_moon(times_Current)
    moonaltazs_Current = moon_Current.transform_to(frame_Current)

    m33altazs_Current = my_Loc_altaz.transform_to(frame_Current)

    #plot it in a pretty graph with a colorbar

    ax.set_title('Astronomy Object Finder')
    ax.plot(delta_midnight, sunaltazs_Current.alt, color='yellow', label='Sun')
    ax.plot(delta_midnight, moonaltazs_Current.alt, color=[0.75]*3, ls='--', label='Moon')
    im = ax.scatter(delta_midnight, m33altazs_Current.alt, c=m33altazs_Current.az, label=obj, lw=0, s=8, cmap='viridis')
    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg, sunaltazs_Current.alt < -0*u.deg, color='0.5', zorder=0)
    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg, sunaltazs_Current.alt < -18*u.deg, color='k', zorder=0)
    
    ax.legend(loc='upper left')
    ax.set_xlim(-12*u.hour, 12*u.hour)
    ax.set_xticks((np.arange(13)*2-12)*u.hour)
    ax.set_ylim(-15*u.deg, 95*u.deg)
    ax.set_xlabel('Hours from ^NOW^ Midnight')
    ax.set_ylabel('Altitude [deg]')
    ax.set_title(obj + ' Transit')
    style = dict(size=10, color='blue')
    ax.text(0,curAlt,curpt,**style)

    plt.tight_layout()


    # Add a colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

#Animate the graph and update every 5 seconds
ani = animation.FuncAnimation(fig, animate, interval=5)
plt.show()
