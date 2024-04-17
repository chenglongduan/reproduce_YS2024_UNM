# 1. STA/LTA PP image: Z_pp_ps
# source wavelet time shift = 0

# raw spectrum
< vz_pp_dt0.rsf sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# raw image
< vz_pp_dt0.rsf sfgrey color=J scalebar=y | sfpen

# high-pass spectrum
< vz_pp_dt0.rsf sfbandpass flo=0.0005 | sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# high-pass image
< vz_pp_dt0.rsf sfbandpass flo=0.0005 | sfgrey color=J scalebar=y | sfpen

# (matlab)horizontal median-filter for surface part
< vz_pp_dt0.rsf sfbandpass flo=0.0005 datapath=./ > vz_pp_dt0_final.rsf
......run matlab......
# check final image
< vz_pp_dt0_final.rsf sfgrey color=J scalebar=y | sfpen

# Kz-X spectrum
< vz_pp_dt0_final.rsf sfspectra all=n | sfgrey color=J max1=0.004 title="PP image" scalebar=y | sfpen
< vz_pp_dt0_final.rsf sfspectra all=n datapath=./ > Kz_vz_pp_dt0_final.rsf
#====================================================================================================

# 2. STA/LTA PS image: R_ps
# source wavelet time shift = 0
# The original mute is better

# raw spectrum
< vr_ps_dt0.rsf sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# raw image
< vr_ps_dt0.rsf sfgrey color=J scalebar=y | sfpen

# high-pass spectrum
< vr_ps_dt0.rsf sfbandpass flo=0.001 | sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# high-pass image
< vr_ps_dt0.rsf sfbandpass flo=0.001 | sfgrey color=J scalebar=y | sfpen

# (matlab)horizontal median-filter for surface part
< vr_ps_dt0.rsf sfbandpass flo=0.001 datapath=./ > vr_ps_dt0_final.rsf
......run matlab......
# check final image
< vr_ps_dt0_final.rsf sfgrey color=J scalebar=y | sfpen
#====================================================================================================

# 3. STA/LTA PS image: T_ps
# source wavelet time shift = 0
# The original mute is better

# raw spectrum
< vt_ps_dt0.rsf sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# raw image
< vt_ps_dt0.rsf sfgrey color=J scalebar=y | sfpen

# high-pass spectrum
< vt_ps_dt0.rsf sfbandpass flo=0.001 | sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# high-pass image
< vt_ps_dt0.rsf sfbandpass flo=0.001 | sfgrey color=J scalebar=y | sfpen

# (matlab)horizontal median-filter for surface part
< vt_ps_dt0.rsf sfbandpass flo=0.001 datapath=./ > vt_ps_dt0_final.rsf
......run matlab......
# check final image
< vt_ps_dt0_final.rsf sfgrey color=J scalebar=y | sfpen
#====================================================================================================

# 4. STA/LTA PS image: R+T

sfmath x=vr_ps_dt0.rsf y=vt_ps_dt0.rsf output="x+y" datapath=./ > R_plus_T.rsf

# raw spectrum
< R_plus_T.rsf sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# raw image
< R_plus_T.rsf sfgrey color=J scalebar=y | sfpen

# high-pass spectrum
< R_plus_T.rsf sfbandpass flo=0.001 | sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" max1=0.006 | sfpen
# high-pass image
< R_plus_T.rsf sfbandpass flo=0.001 | sfgrey color=J scalebar=y | sfpen

# (matlab)horizontal median-filter for surface part
< R_plus_T.rsf sfbandpass flo=0.001 datapath=./ > R_plus_T_final.rsf
......run matlab......
# check final image
< R_plus_T_final.rsf sfgrey color=J scalebar=y | sfpen

# Kz-X spectrum
< R_plus_T_final.rsf sfspectra all=n | sfgrey color=J max1=0.004 title="PS image" scalebar=y | sfpen
< R_plus_T_final.rsf sfspectra all=n datapath=./ > Kz_R_plus_T_final.rsf

