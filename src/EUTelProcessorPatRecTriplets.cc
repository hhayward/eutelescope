#include "EUTelProcessorPatRecTriplets.h"
/**  EUTelProcessorPatRecTriplets
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of hits in local coordinate frame of telescope planes. 
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of track candidates.
 */
//If you are using static variables then for the most part you are doing something wrong. 
//TO DO: Create new Histograming procedure! 

using namespace eutelescope;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorPatRecTriplets::_histName::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelProcessorPatRecTriplets::_histName::_numberOfHitOnTrackCandidateHistName = "NumberOfHitsOnTrackCandidate";
std::string EUTelProcessorPatRecTriplets::_histName::_HitOnTrackCandidateHistName = "HitsOnTrackCandidate";
std::string EUTelProcessorPatRecTriplets::_histName::_chi2CandidateHistName = "chi2CandidateHistName";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorPatRecTriplets::EUTelProcessorPatRecTriplets() :
Processor("EUTelProcessorPatRecTriplets"),
_aidaHistoMap1D(),
_hitInputCollectionName("HitCollection"),
_trackCandidateHitsOutputCollectionName("TrackCandidatesCollection"),
_histoInfoFileName("histoinfo.xml"),
_trackFitter(0),
_doubletDistCut(0,0),
_tripletSlopeCuts(0,0),
_nProcessedRuns(0),
_nProcessedEvents(0),
_eBeam(-1.),
_qBeam(-1.)
{
	//The standard description that comes with every processor 
	_description = "EUTelProcessorPatRecTriplets preforms track pattern recognition.";

	// TrackerHit input collection
	registerInputCollection(LCIO::TRACKERHIT,"HitInputCollectionName","Input hits collection name",_hitInputCollectionName,std::string("HitCollection"));

	// Track candidate hits output collection
	registerOutputCollection(LCIO::TRACK,"TrackCandHitOutputCollectionName","Output track candidates hits collection name",_trackCandidateHitsOutputCollectionName,std::string("TrackCandidateHitCollection"));

	registerOptionalParameter("DoubletDistCut", "Doublet distance cuts", _doubletDistCut, FloatVec()); 

	//This is a another cut we do after we find the suspected track from a single particle
	//If this number of hits are shared with another possible track then we discard one of the tracks. 
	registerOptionalParameter("TripletSlopeCuts", "Triplet slope difference which is allowed to create track ",
	_tripletSlopeCuts, FloatVec());

	//This cut is used during the search for hits. If we propagate to a plane and the nearest hit to the point of intersection is greater than this number then we assume the hit can not come from this track. 
	registerOptionalParameter("DoubletCenDistCut","Doublet hit acceptance distance for central plane ", _doubletCenDistCut,FloatVec() );
	registerOptionalParameter("TripletConnectDistCut","The distance cut to allow the triplets to be associated with each other", _tripletConnectDistCut,FloatVec() );

	//This is need if we have a magnetic field to determine curvature
	registerOptionalParameter("BeamEnergy", "Beam energy [GeV]", _eBeam, static_cast<double> (4.0));

	//This is needed if we have a magnetic field to determine curvature. 
	registerOptionalParameter("BeamCharge", "Beam charge [e]", _qBeam, static_cast<double> (-1));

	// Histogram information. This is the xml file that specifies the name and bin size etc of the create histograms.
	registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

	//This is the planes at which we will begin to look for tracks from. We start from this plane and move forward looking to see if hits can be added to the possible track. 
	registerOptionalParameter("PlanesToCreateSeedsFrom", "This is the planes you want to create seeds from", _createSeedsFromPlanes,IntVec());

	//This is planes that we should not look for hits or create a state. Effectively this removes this plane from the analysis. However scattering due the material is still taken into account
	registerOptionalParameter("ExcludePlanes", "This is the planes that will not be included in analysis", _excludePlanes ,IntVec());

	//This initial distance to the first plane is needed if we are in a magnetic field. Since the particle will begin to curve before we reach the first plane. Therefore to accurately model the particles trajectory we need to take this into account. 
	registerOptionalParameter("InitialDisplacement", "This is the initial distance the particle must travel to reach the first plane", _initialDisplacement ,float(0));

	//This specifies if the planes are strip or pixel sensors.
	registerOptionalParameter("planeDimensions", "This is a number 1(strip sensor) or 2(pixel sensor) to identify the type of detector. Must be in z order and include all planes.", _planeDimension, IntVec());

}
//This is the inital function that Marlin will run only once when we run jobsub
void EUTelProcessorPatRecTriplets::init(){

	try{
		_nProcessedRuns = 0;
		_nProcessedEvents = 0;
		std::string name = EUTELESCOPE::GEOFILENAME;
		geo::gGeometry().initializeTGeoDescription(name,false);
		
		geo::gGeometry().initialisePlanesToExcluded(_excludePlanes);//We specify the excluded planes here since this is rather generic and can be used by other processors
		geo::gGeometry().setInitialDisplacementToFirstPlane(_initialDisplacement);//We specify this here so we can access it throughout this processor. 
		
		streamlog_out(MESSAGE5) << "These are the planes you will create a state from. Mass inbetween states will be turned to scatterers in GBLTrackProcessor." << std::endl;
		for(size_t i =0 ; i < geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size(); ++i)
		{
			streamlog_out(MESSAGE5) << geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) << "  ";
		}
		streamlog_out(MESSAGE5) << std::endl;

		_trackFitter = new EUTelPatRecTriplets();//This is the class that contains all the functions that do the actual work
		
		_trackFitter->setDoubletDistCut(_doubletDistCut);
		_trackFitter->setTripletSlopeCuts(_tripletSlopeCuts);
		_trackFitter->setDoubletCenDistCut(_doubletCenDistCut);
        _trackFitter->setTripletConnectDistCut(_tripletConnectDistCut);
		_trackFitter->setPlanesToCreateSeedsFrom(_createSeedsFromPlanes);
		_trackFitter->setBeamMomentum(_eBeam);
		_trackFitter->setBeamCharge(_qBeam);
		_trackFitter->setAutoPlanestoCreateSeedsFrom();//If the user has not specified which planes to seed from the the first plane is used
		_trackFitter->testUserInput();//Here we check that the user has provided the correct data. This is the most likey place to throw and exception.
		bookHistograms();		// Book histograms. Yet again this should be replaced. TO DO:Create better histogram method.
	}
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;

	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in init function of pattern recognition" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorPatRecTriplets::processRunHeader(LCRunHeader* run) {

	std::auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(run));
	header->addProcessor(type());//Add what processor has acted to collection here. 

	// this is the right place also to check the geometry ID. This is a
	// unique number identifying each different geometry used at the
	// beam test. The same number should be saved in the run header and
	// in the xml file. If the numbers are different, warn the user.

	if (header->getGeoID() == 0)
		streamlog_out(WARNING0) << "The geometry ID in the run header is set to zero." << std::endl
		<< "This may mean that the GeoID parameter was not set" << std::endl;


	if ((unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID()) {
		streamlog_out(WARNING5) << "Error during the geometry consistency check: " << std::endl
				<< "The run header says the GeoID is " << header->getGeoID() << std::endl
				<< "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
	}
		_nProcessedRuns++;
}

void EUTelProcessorPatRecTriplets::processEvent(LCEvent* evt)
{
	UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

	try{
		EUTelEventImpl* event = static_cast<EUTelEventImpl*> (evt); //Change the LCIO object to EUTel object. This is a simple way to extend functionality of the object.
		_trackFitter->setEventNumber(_nProcessedEvents);//This is so we can use the event number with this class. 
		// Do not process last event. For unknown events just a warning will do. 
		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		} else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber()<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
		LCCollection* hitMeasuredCollection = NULL;
		hitMeasuredCollection = evt->getCollection(_hitInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " <<_hitInputCollectionName << " retrieved" << std::endl;

		// this will only be entered if the collection is not available
		if ( hitMeasuredCollection == NULL) {
			streamlog_out(DEBUG2) << "EUTelProcessorPatRecTriplets :: processEvent() hitMeasuredCollection is void" << std::endl;
			throw marlin::SkipEventException(this);
		}

		// Prepare hits for track finder
		EVENT::TrackerHitVec allHitsVec;
		_trackFitter->clearFinalTracks(); //This is to clear the vector of tracks from the last event.

		for(int iHit = 0; iHit < hitMeasuredCollection->getNumberOfElements(); iHit++) 
		{
			TrackerHitImpl* hit = static_cast<TrackerHitImpl*>(hitMeasuredCollection->getElementAt(iHit));
                	int sensorID = hitDecoder(static_cast<IMPL::TrackerHitImpl*>(hit))["sensorID"];
			if ( sensorID >= 0 ) allHitsVec.push_back(hit);
		}

		//TODO: this should never happen
		if(allHitsVec.empty()) throw lcio::Exception("No hits!");

		_trackFitter->setHitsVec(allHitsVec);//Will set data member of class _trackFitter. TO DO: We could do this within findHitsOrdervec.  

//		_trackFitter->printHits();//We print the hits to screen for debugging
//		streamlog_out(DEBUG2) << "Now map hits to planes" << std::endl;

		_trackFitter->setHitsVecPerPlane();//Create map Sensor ID(non excluded)->HitsVec (Using geometry)

//		streamlog_out(DEBUG2) << "Test hits on the planes" << std::endl;

		_trackFitter->testHitsVecPerPlane();//tests the size of the map and does it contain hits
		streamlog_out( DEBUG1 ) << "Trying to find tracks..." << std::endl;
		_trackFitter->createTriplets();

		_trackFitter->getTracks();
//		plotHistos(tracks);

//		outputLCIO(evt,tracks);

		_nProcessedEvents++;
	}
	catch (DataNotAvailableException e) {
		streamlog_out(WARNING2) << " Collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in process function of pattern recognition" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorPatRecTriplets::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>& tracks)
{
    if(!tracks.empty()){
        for(unsigned int i=0 ; i< tracks.size(); i++){
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" track for event " << evt->getEventNumber() <<".  Track number  " << i <<std::endl;
            tracks.at(i).print();
        }
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        reader.getColVec(tracks, evt, _trackCandidateHitsOutputCollectionName);
    }
}

//TO DO: find a more generic way to plot histograms
void EUTelProcessorPatRecTriplets::plotHistos( std::vector<EUTelTrack>& trackCandidates)  {

	const int nTracks = trackCandidates.size( );
 
	static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> fill( nTracks );
	streamlog_out( MESSAGE2 ) << "Event #" << _nProcessedEvents << std::endl;
	int numberOfHits =0;
	for (size_t i = 0; i< trackCandidates.size( ) ; ++i ) {//loop over all tracks
		for(size_t j = 0; j <trackCandidates[i].getStates().size(); ++j){//loop over all states 
			if(!trackCandidates[i].getStates()[j].getStateHasHit()){//We only ever store on hit per state
				continue;
			}
			numberOfHits++;
			int sensorID = static_cast<int>(trackCandidates[i].getStates().at(j).getLocation());//since we store sensor ID in Z0
			static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_HitOnTrackCandidateHistName ] ) -> fill( sensorID );
		}
	streamlog_out( MESSAGE1 ) << "Track hits end:==============" << std::endl;
	static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> fill(numberOfHits);
		}
}

void EUTelProcessorPatRecTriplets::check(LCEvent * /*evt*/) {
    // nothing to check here
}

void EUTelProcessorPatRecTriplets::end() {
    streamlog_out(MESSAGE9) <<"The average number of tracks per event: " << static_cast<float>(_trackFitter->getNumberOfTracksAfterPruneCut())/static_cast<float>(_nProcessedEvents) <<std::endl; 
    delete _trackFitter;
    streamlog_out(MESSAGE9) << "EUTelProcessorPatRecTriplets::end()  " << name()
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << " av.tracks : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberTracksCandidatesHistName ] ) -> mean()
            << " track candidates : " << static_cast < AIDA::IHistogram1D* > ( _aidaHistoMap1D[ _histName::_numberOfHitOnTrackCandidateHistName ] ) -> allEntries()
            << std::endl;
}
//TO DO: Create a better way of booking histograms.
void EUTelProcessorPatRecTriplets::bookHistograms() {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out(DEBUG) << "Booking histograms..." << std::endl;

        std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
        EUTelHistogramInfo    * histoInfo;
        bool                    isHistoManagerAvailable;

        try {
            isHistoManagerAvailable = histoMgr->init( );
        } catch ( std::ios::failure& e ) {
            streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
                    << "Continuing without histogram manager using default settings"    << std::endl;
            isHistoManagerAvailable = false;
        } catch ( marlin::ParseException& e ) {
            streamlog_out( ERROR5 ) << e.what( ) << "\n"
                    << "Continuing without histogram manager using default settings" << std::endl;
            isHistoManagerAvailable = false;
        }
        
        const int tracksNBin = 20;    
        const double tracksXMin = -0.5;
        const double tracksXMax = 19.5;

        histoInfo = histoMgr->getHistogramInfo(_histName::_numberTracksCandidatesHistName);        
        int NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : tracksNBin;
        double XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : tracksXMin;
        double XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : tracksXMax;

        // Number of track candidates from pattern recognition step
        AIDA::IHistogram1D * numberTracksCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_numberTracksCandidatesHistName, NBin,
                XMin, XMax);

        if (numberTracksCandidates) {
            numberTracksCandidates->setTitle("Number of track candidates;N tracks;N events");
            _aidaHistoMap1D.insert(make_pair(_histName::_numberTracksCandidatesHistName, numberTracksCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberTracksCandidatesHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
	        }
			const int chiNbins = 5000;
			const double chiXmin = 0;
			const double chiXmax = 500000;												


				///////////////////////////////////////////////////////////////////////////////////////////////////////////////Chi2 create plot. Useful to determine the behaviour of the pattern recognition
        histoInfo = histoMgr->getHistogramInfo(_histName::_chi2CandidateHistName);
//      NBin =  histoInfo->_xBin; 
//      XMin =  histoInfo->_xMin; 
//      XMax = histoInfo->_xMax;
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : chiNbins;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :chiXmin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : chiXmax;

        AIDA::IHistogram1D * chi2 =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_chi2CandidateHistName, NBin,
                XMin, XMax);

            numberTracksCandidates->setTitle("Number of track candidates;N tracks;N events");
            _aidaHistoMap1D.insert(make_pair(_histName::_chi2CandidateHistName, chi2));
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				  
				
         
        const int hitsNBin = 10;
        const double hitsMin = -0.5;
        const double hitsMax = 9.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_numberOfHitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : hitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : hitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : hitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * numberHitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_numberOfHitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (numberHitsOnTrackCandidates) {
            numberHitsOnTrackCandidates->setTitle("Number of hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_numberOfHitOnTrackCandidateHistName, numberHitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_numberOfHitOnTrackCandidateHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

        const int    nhitsNBin = 40;
        const double nhitsMin = -0.5;
        const double nhitsMax =39.5;
        
        histoInfo = histoMgr->getHistogramInfo(_histName::_HitOnTrackCandidateHistName);        
        NBin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : nhitsNBin;
        XMin = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin : nhitsMin;
        XMax = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : nhitsMax;
        
        // Number of hit per track candidate
        AIDA::IHistogram1D * HitsOnTrackCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_histName::_HitOnTrackCandidateHistName, NBin,
                XMin, XMax);

        if (HitsOnTrackCandidates) {
            HitsOnTrackCandidates->setTitle("hits on track candidates;N hits;N tracks");
            _aidaHistoMap1D.insert(make_pair(_histName::_HitOnTrackCandidateHistName, HitsOnTrackCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_histName::_HitOnTrackCandidateHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }
       
    } catch (lcio::Exception& e) {
        streamlog_out(WARNING2) << "Can't allocate histograms. Continue without histogramming" << std::endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

