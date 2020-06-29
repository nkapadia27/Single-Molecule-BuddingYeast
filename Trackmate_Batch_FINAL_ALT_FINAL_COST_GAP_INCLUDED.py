from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from ij import IJ
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.features.FeatureAnalyzer as FeatureAnalyzer
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory as SpotContrastAndSNRAnalyzerFactory
import fiji.plugin.trackmate.action.ExportStatsToIJAction as ExportStatsToIJAction
import fiji.plugin.trackmate.io.TmXmlReader as TmXmlReader
import fiji.plugin.trackmate.action.ExportTracksToXML as ExportTracksToXML
import fiji.plugin.trackmate.io.TmXmlWriter as TmXmlWriter
import fiji.plugin.trackmate.features.ModelFeatureUpdater as ModelFeatureUpdater
import fiji.plugin.trackmate.features.SpotFeatureCalculator as SpotFeatureCalculator
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzer as SpotContrastAndSNRAnalyzer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory as SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.features.spot.SpotRadiusEstimatorFactory as SpotRadiusEstimatorFactory
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer
import fiji.plugin.trackmate.features.track.TrackSpotQualityFeatureAnalyzer as TrackSpotQualityFeatureAnalyzer
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.features.track.TrackLocationAnalyzer as TrackLocationAnalyzer
import fiji.plugin.trackmate.features.track.TrackBranchingAnalyzer as TrackBranchingAnalyzer
import fiji.plugin.trackmate.features.track.TrackIndexAnalyzer as TrackIndexAnalyzer
import fiji.plugin.trackmate.SpotCollection as SpotCollection
import fiji.plugin.trackmate.util.TMUtils as TMUtils
import csv
import os
import math
import itertools as itto



# Get currently selected image
#imp = WindowManager.getCurrentImage()
dir_tiff = "C:/Users/Reyes-Lamothe Lab/Pictures/Microscopy Data/Nitin/PALM Yeast/20190709/ZEY127_HALO_8s/Tiff Files/Timelapses"
link_dist = 3.0 
gap_link_dist = 5.0 
split_dist = 5.0 
merg_dist = 5.0 
int_thresh = 500.
cost_q = 0.3
for (root, dirs, files) in os.walk(dir_tiff):
    files = [ fi for fi in files if fi.endswith(".tif") ]
    
    num_files = len(files)
    print(files)
    print (root)
save_analysis_dir = dir_tiff + '/' + 'Analysis_COST' + '_' + 'thr' + str(int_thresh) + '_' + 'lnk' + str(link_dist) + '_' + 'gp' + str(gap_link_dist) + '_' + 'spl' + str(split_dist) + '_' + 'mrg' + str(merg_dist) + '_' + 'cost_INT' + str(cost_q)
os.mkdir(save_analysis_dir)
for z in range(num_files):
    str_name_fil = root + '/' + files [z]
    print(str_name_fil)
    imp = IJ.openImage(str_name_fil)
    #imp.show()
    dims = imp.getDimensions()
    imp.setDimensions( dims[ 2 ], dims[ 4 ], dims[ 3 ] )

    #----------------------------
    # Create the model object now
    #----------------------------

    # Some of the parameters we configure below need to have
    # a reference to the model at creation. So we create an
    # empty model now.

    model = Model()

    # Send all messages to ImageJ log window.
    model.setLogger(Logger.IJ_LOGGER)



    #------------------------
    # Prepare settings object
    #------------------------

    settings = Settings()
    settings.setFrom(imp)

    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = LogDetectorFactory()
    settings.detectorSettings = {
        'DO_SUBPIXEL_LOCALIZATION' : True,
        'RADIUS' : 1.25,
        'TARGET_CHANNEL' : 1,
        'THRESHOLD' : int_thresh,
        'DO_MEDIAN_FILTERING' : True,
    }
    settings.dx = 100.0
    settings.dy = 100.0


    # Configure spot filters - Classical filter on quality
    #filter1 = FeatureFilter('QUALITY', 30, True)
    #settings.addSpotFilter(filter1)

    # Configure tracker - We want to allow merges and fusions
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap() # almost good enough
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
    settings.trackerSettings['MAX_FRAME_GAP'] = 1
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = link_dist
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = gap_link_dist
    settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = split_dist
    settings.trackerSettings['MERGING_MAX_DISTANCE'] = merg_dist
    settings.trackerSettings['LINKING_FEATURE_PENALTIES'] = {'QUALITY' : cost_q}
    settings.trackerSettings['GAP_CLOSING_FEATURE_PENALTIES'] = {'QUALITY' : cost_q} 
    # Configure track analyzers - Later on we want to filter out tracks
    # based on their displacement, so we need to state that we want
    # track displacement to be calculated. By default, out of the GUI,
    # not features are calculated.

    # The displacement feature is provided by the TrackDurationAnalyzer.



    # Configure track filters - We want to get rid of the two immobile spots at
    # the bottom right of the image. Track displacement must be above 10 pixels.

    #filter2 = FeatureFilter('TRACK_DISPLACEMENT', 10, True)
    #settings.addTrackFilter(filter2)
    settings.addSpotAnalyzerFactory(SpotRadiusEstimatorFactory())
    settings.addSpotAnalyzerFactory(SpotIntensityAnalyzerFactory())
    settings.addSpotAnalyzerFactory(SpotContrastAndSNRAnalyzerFactory())
	
    
    # Add an analyzer for some track features, such as the track mean speed.
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.addTrackAnalyzer(TrackSpeedStatisticsAnalyzer())
    settings.addTrackAnalyzer(TrackSpotQualityFeatureAnalyzer())
    settings.addTrackAnalyzer(TrackLocationAnalyzer())
    settings.addTrackAnalyzer(TrackBranchingAnalyzer())
    settings.addTrackAnalyzer(TrackIndexAnalyzer())
    #settings.initialSpotFilterValue = 1

    print(str(settings))

    #-------------------
    # Instantiate plugin
    #-------------------

    trackmate = TrackMate(model, settings)
    trackmate.setNumThreads(8)
    #--------
    # Process
    #--------

    ok = trackmate.checkInput()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))

    ok = trackmate.process()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))


    #----------------
    # Display results
    #----------------

    #model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')

    selectionModel = SelectionModel(model)
    #displayer =  HyperStackDisplayer(model, selectionModel, imp)
    #displayer.render()
    #displayer.refresh()

    # The feature model, that stores edge and track features.
    fm = model.getFeatureModel()
    #print (len(model.getTrackModel().trackIDs(True)))
    Track_IDs = [None]*0
    mean_sp = [None]*0
    med_sp = [None]*0
    min_sp = [None]*0
    max_sp = [None]*0
    std_sp = [None]*0
    mean_q = [None]*0
    med_q_tr = [None]*0
    min_q_tr = [None]*0
    max_q_tr = [None]*0
    std_q_tr = [None]*0
    x_lc = [None]*0
    y_lc = [None]*0
    mn_int = [None]*0
    inten = [None]*0
    tr_dur = [None]*0
    tr_start = [None]*0
    tr_fin = [None]*0
    spt_tr = [None]*0
    spt_widt = [None]*0
    tr_charact = [None]*0
    x_tr = [None]*0
    y_tr = [None]*0
    spt_all_x = [None]*0
    spt_all_y = [None]*0
    tr_identifi = [None]*0
    tr_fram = [None]*0
    spt_m = model.getSpots()
    #tot_spts = spt_m.getNSpots(True)
    for id in model.getTrackModel().trackIDs(True):

        # Fetch the track feature from the feature model.
        v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
        med_v = fm.getTrackFeature(id, 'TRACK_MEDIAN_SPEED')
        min_v = fm.getTrackFeature(id, 'TRACK_MIN_SPEED')
        max_v = fm.getTrackFeature(id, 'TRACK_MAX_SPEED')
        std_v = fm.getTrackFeature(id, 'TRACK_STD_SPEED')
        q = fm.getTrackFeature(id, 'TRACK_MEAN_QUALITY')
        med_q = fm.getTrackFeature(id, 'TRACK_MEDIAN_QUALITY')
        min_q = fm.getTrackFeature(id, 'TRACK_MIN_QUALITY')
        max_q = fm.getTrackFeature(id, 'TRACK_MAX_QUALITY')
        std_q = fm.getTrackFeature(id, 'TRACK_STD_QUALITY')
        dura = fm.getTrackFeature(id, 'TRACK_DURATION')
        start_tr = fm.getTrackFeature(id, 'TRACK_START')
        fin_tr = fm.getTrackFeature(id, 'TRACK_STOP')
        spts = fm.getTrackFeature(id, 'NUMBER_SPOTS')
        tr_identif = fm.getTrackFeature(id,'TRACK_ID')
        #print(spts)
        identi = [int(id)]*int(spts)
        
        #print(identi)
        #x_loc = fm.getTrackFeature(id, 'X_LOCATION')
        #print(x_loc)
        #print(x_loc)
        #y_loc = fm.getTrackFeature(id, 'Y_LOCATION')
        #model.getLogger().log('')
        #model.getLogger().log('Track ' + str(id) + ': mean velocity = ' + str(v) + ' ' + model.getSpaceUnits() + '/' + model.getTimeUnits())

        track = model.getTrackModel().trackSpots(id)
        track_int = [None]*0
        track_x = [None]*0
        track_y = [None]*0
        track_coord = [None]*0
        inten2 = [None]*0
        track_shp = [None]*0
        fram = [None]*0
        #identi = [None]*0
        for spot in track:
            sid = spot.ID()
            # Fetch spot features directly from spot.
            x=spot.getFeature('POSITION_X')
            y=spot.getFeature('POSITION_Y')
            t=spot.getFeature('FRAME')
            #q=spot.getFeature('QUALITY')
            snr=spot.getFeature('SNR')
            tot_int_spot=spot.getFeature('TOTAL_INTENSITY')
            wid=spot.getFeature('ESTIMATED_DIAMETER')
            track_int.append(tot_int_spot)
            #inten2.append(mean_int_spot)
            track_x.append(x)
            track_y.append(y)
            fram.append(t)
            #track_coord.append(t, x, y)
            track_shp.append(wid)
        #print(len(track_x[0]))
    
        #dattr = track_x
        #print(indenti_2)
        Track_IDs.append(id)# = id
        mean_sp.append(v)
        med_sp.append(med_v)# =
        min_sp.append(min_v) #= min_v
        max_sp.append(max_v)
        std_sp.append(std_v)#= max_v
        mean_q.append(q)
        med_q_tr.append(med_q)#= q
        min_q_tr.append(min_q)
        max_q_tr.append(max_q)
        std_q_tr.append(std_q)#= min_q
        x_loc = sum(track_x)/len(track_x)
        y_loc = sum(track_y)/len(track_y)
        x_lc.append(x_loc) #= x_loc
        y_lc.append(y_loc) #= y_loc
        mean_track_int = sum(track_int)/len(track_int)
        mn_int.append(mean_track_int) #= mean_int
        inten = inten + track_int
        tr_dur.append(dura)
        tr_start.append(start_tr)
        tr_fin.append(fin_tr)
        spt_tr.append(spts)
        spt_widt.append(sum(track_shp)/len(track_shp))
        tr_fram = tr_fram + fram
        x_tr = x_tr + track_x
        y_tr = y_tr + track_y
        tr_identifi = tr_identifi + identi
   	spt_all_fin = [tr_identifi, tr_fram, x_tr, y_tr, inten]
	spt_all_fin_2 = [[row[i] for row in spt_all_fin]
                             for i in range(len(spt_all_fin[0]))]
         
    data_final = [Track_IDs, spt_tr, spt_widt, mean_sp, max_sp, min_sp, med_sp, std_sp, mean_q, max_q_tr, min_q_tr, med_q_tr, std_q_tr, tr_dur, tr_start, tr_fin, x_lc, y_lc]
    data_final_2 = [[row[i] for row in data_final]
                             for i in range(len(data_final[0]))]
	
	
	
	

    
    dir_name_save = save_analysis_dir + '/' + files[z] + '_' + 'tracks'
    dir_name_save_inten = save_analysis_dir + '/' + files[z] + '_' + 'spots'
    str_save = dir_name_save + '.csv'
    str_save_inten = dir_name_save_inten + '.csv'
    with open(str_save, "w") as f:
        writer = csv.writer(f,delimiter=',')
        writer.writerows(data_final_2)
    with open(str_save_inten, "w") as f:
        writer = csv.writer(f,delimiter=',')
        writer.writerows(spt_all_fin_2)
    IJ.log('Success')
    del data_final_2
    del spt_all_fin_2
    imp.close()
    imp.flush()

