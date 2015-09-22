/* To create a new analysis follow the series of steps below:
 * 1) Create a new set of histograms in this processor and link the histogram's pointers to the state's locations; i.e the plane number.
 * 		Histograms are made using the input to the processor and then matched using the input from the states location.
 * 2) Create a new input member function to set this in the class (EUTelTrackAnalysis).
 * 3) Now pass the track from this processor to EUTelTrackAnalysis via a function as shown in processEvent below.
 * 4)You now have the trackand histogram. Do the analysis and output to that histogram or anyone oyu want.   */
#include "EUTelProcessorTrackAnalysis.h"
#include <AIDA/IAxis.h>      

// eutelescope geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

using namespace eutelescope;

EUTelProcessorTrackAnalysis::EUTelProcessorTrackAnalysis() :
Processor("EUTelProcessorTrackAnalysis"){
	
	registerInputCollection(LCIO::TRACK, "TrackInputCollectionName", "Input track collection name",_trackInputCollectionName,std::string("TrackCandidatesCollection"));
	registerOptionalParameter("SensorIDs", "A vector of the sensor IDs to plot",_sensorIDs,IntVec());
    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

	
	
}


void EUTelProcessorTrackAnalysis::init(){

	//Helen trying to access pitch information
 // std::cout<<"Helen: lets try to retrieve the geo information"     <<std::endl;
 // 	for( EVENT::IntVec::iterator it = _sensorIDs.begin(); it != _sensorIDs.end(); ++it ) 
 // 	  {
 // 	    std::cout<<"Helen: sensor id = "<<*it<<std::endl;
 // 	    try
 // 	      {
	
 // 	      }
 // 	catch(std::runtime_error& e)
 // 		{
 // 			streamlog_out ( ERROR0 ) << "hELEN could not retrieve plane " << *it << std::endl;
 // 			streamlog_out ( ERROR0 ) << e.what() << std::endl;
 // 			//throw StopProcessingException(this);
 // 		}
 // 	  }
	try{
		initialiseResidualVsPositionHistograms();
		initialiseHitMapHistograms();
		initialiseEfficiencyVsPositionHistograms();
		initialiseGeoInfo();//fill array of pitch info here
		//Some initialised in the constructor in part 2.
		EUTelTrackAnalysis*	analysis = new EUTelTrackAnalysis(_mapFromSensorIDToHistogramX,_mapFromSensorIDToHistogramY,_mapFromSensorHitMap,_mapFromSensorIDToEfficiencyX,_mapFromSensorIDToEfficiencyY,_mapFromSensorIDToGloIncXZ,_mapFromSensorIDToGloIncYZ,_beamEnergy, _mapFromSensorPitchX,_mapFromSensorPitchY); 

		//Others here.
		analysis->setSensorIDTo2DPValuesWithPosition(_mapFromSensorIDToPValueHisto);
		analysis->setSensorIDToPValuesVsIncidenceAngleYZ(_mapFromSensorIDToPValuesVsIncidenceYZ);
		analysis->setSensorIDToPValuesVsIncidenceAngleXZ(_mapFromSensorIDToPValuesVsIncidenceXZ);

		analysis->setPValueBeamEnergy(_pValueVsBeamEnergy);
		_analysis = analysis;
	}catch(...){	
		streamlog_out(MESSAGE9)<<"There is an unknown error in EUTelProcessorTrackAnalysis-init()" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}


}

void EUTelProcessorTrackAnalysis::processEvent(LCEvent * evt){
	try{
		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
        streamlog_out(DEBUG2) << "Collection contains data! Continue!" << std::endl;
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO(); streamlog_out(DEBUG2) << "Collection contains data! Continue! line 53" << std::endl;
        streamlog_out(DEBUG2) << "_trackInputCollectionName = " <<_trackInputCollectionName<<std::endl;
        std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackInputCollectionName);
        streamlog_out(DEBUG2) << "Collection contains data! Continue! line 54: tracks.size() = " << tracks.size()<<std::endl;
        for (int iTrack = 0; iTrack < tracks.size(); ++iTrack){
   //         track.print();
            EUTelTrack track = tracks.at(iTrack); 
	    //std::cout<<"--Helen: track "<<iTrack<<" out of "<<tracks.size()<<std::endl;
            _analysis->plotResidualVsPosition(track);
	    _analysis->plotHitMap(track);
            _analysis->plotEfficiencyVsPosition(track,_sensorIDs);	
            _analysis->plotIncidenceAngles(track);
            if(track.getChi2()/track.getNdf() < 5.0){
                _analysis->plotBeamEnergy(track);
                _analysis->plotPValueVsBeamEnergy(track);
            }//if(track.getChi2()/track.getNdf() < 5.0){
            _analysis->plotPValueWithPosition(track);
            _analysis->plotPValueWithIncidenceAngles(track);
           _analysis->setTotNum(track);

        }//for (int iTrack = 0; iTrack < tracks.size(); ++iTrack){
        }catch (DataNotAvailableException e) {
//		streamlog_out(WARNING2) << " Collection not available" << std::endl;
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
		streamlog_out(MESSAGE9)<<"Unknown exception in process function of track analysis" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}

	
}

void EUTelProcessorTrackAnalysis::end(){

  streamlog_out(DEBUG2) <<" HELEN here"<<std::endl;
  
  //test i can add stuff here
  //1D histogram of efficiency value - one per sensor
  std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
  bool isHistoManagerAvailable;
  std::stringstream sstm;
  std::string elementEffDistHistName;
  std::string elementEffDistHistLooseName;
  std::string histTitle;
  // isHistoManagerAvailable = histoMgr->init( );
  //streamlog_out(DEBUG2) <<" HELEN inside postprocessing, isHistoManagerAvailable = "<< isHistoManagerAvailable<<std::endl;
  for(size_t i = 0; i < _sensorIDs.size() ; ++i){
    sstm << "elementEffDist" << _sensorIDs.at(i);
    elementEffDistHistName = sstm.str();
    sstm.str(std::string());
    sstm << "elementEffDist " <<  _sensorIDs.at(i);
    histTitle = sstm.str();
    sstm.str(std::string(""));
    histoInfo = histoMgr->getHistogramInfo(elementEffDistHistName);
    AIDA::IHistogram1D * elementEffDist = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(elementEffDistHistName, 101, 0.0, 1.01); 
    streamlog_out(DEBUG2) <<" HELEN _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->allEntries(); = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->allEntries()<<std::endl;

    //loose acceptance
    sstm << "elementEffDistLoose" << _sensorIDs.at(i);
    elementEffDistHistLooseName = sstm.str();
    sstm.str(std::string());
    sstm << "elementEffDistLoose " <<  _sensorIDs.at(i);
    histTitle = sstm.str();
    sstm.str(std::string(""));
    histoInfo = histoMgr->getHistogramInfo(elementEffDistHistLooseName);
    AIDA::IHistogram1D * elementEffDistLoose = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(elementEffDistHistLooseName, 101, 0.0, 1.01); 
    // streamlog_out(DEBUG2) <<" HELEN _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->allEntries(); = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->allEntries()<<std::endl;
    if(_sensorIDs.at(i)>10){
    //loop over x and y
    for(int x = 0; x < _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->xAxis().bins(); ++x){
      if (x>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(-3)
	  && x<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(7)){
	    for(int y = _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(-4.5); y < _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->yAxis().bins(); ++y){
	      if (y>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(4.5))continue;
	      if( _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->binEntries(x,y)>0){
		//if x is between -3 and 7 and y between +/-4.5
		
		if (x>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(-2)
		    && x<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(6)
		    && y>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(-3.5)
		    && y<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(3.5)){
		  elementEffDist->fill(_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->binHeight(x,y));
		  if(_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->binHeight(x,y)<0.8){
		    //	    std::cout<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->binHeight(x,y)<<" _mapFromSensorIDToEfficiencyX["<<_sensorIDs.at(i)<<"]->coordToIndexX(-2) = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(-2)<<" _mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(6) = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(6)<<" x = "<<x<<std::endl;
		    //	    std::cout<<"_mapFromSensorIDToEfficiencyX["<<_sensorIDs.at(i)<<"]->coordToIndexY(-3.5) = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(-3.5)<<" _mapFromSensorIDToEfficiencyX["<<_sensorIDs.at(i)<<"]->coordToIndexY(3.5) = "<<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(3.5)<<" y = "<<y<<std::endl;
		  }
		}
		
		if (x>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(-3)
		    && x<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexX(7)
		    && y>_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(-4.5)
		    && y<_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->coordToIndexY(4.5)){
		  elementEffDistLoose->fill(_mapFromSensorIDToEfficiencyX[_sensorIDs.at(i)]->binHeight(x,y));
		  
		}
		
	      }
	//streamlog_out(DEBUG2) <<" HELEN still here"<<std::endl;
	    }
	  }
	  }
    //if entry
    //fill elementEffDist

    std::cout<<"Tight acceptance cut [-2<x<6, -3.5<y<3.5] Sensor "<< _sensorIDs.at(i) <<"   Efficiency = "<<elementEffDist->mean()<<" rms = "<<elementEffDist->rms()<<std::endl;
    std::cout<<"Loose acceptance cut [-3<x<7, -4.5<y<4.5] Sensor "<< _sensorIDs.at(i) <<"   Efficiency = "<<elementEffDistLoose->mean()<<" rms = "<<elementEffDistLoose->rms()<<std::endl;
    }
    }
  
    ///Will print the final results of the analysis
    _analysis->print();

}
void	EUTelProcessorTrackAnalysis::initialiseEfficiencyVsPositionHistograms(){
	int NBinX;
	double MinX;
	double MaxX;
	int NBinY;
	double MinY;
	double MaxY;
	double MinZ;
	double MaxZ;

	std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
    bool isHistoManagerAvailable;
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
	isHistoManagerAvailable = false;


	std::stringstream sstm;
	std::string effGblFitHistName;
	std::string histTitle;


	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
		/////////////////////////////////////////////////////////////////////////////XY efficiency plots with position
		sstm << "EfficienyX" << _sensorIDs.at(i);
		effGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "EfficienysX. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(effGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 200;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 200;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 10;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitX =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(effGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitX) {
				residGblFitX->setTitle(histTitle);
				_mapFromSensorIDToEfficiencyX.insert(std::make_pair(_sensorIDs.at(i), residGblFitX));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (effGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "EfficienyY" << _sensorIDs.at(i);
		effGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "EfficienysY. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(effGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 50;//every 500 micron there is a bin
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.0 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.25;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 50;//every 500 micron there is a bin
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -0.0;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 0.05;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitY =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(effGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitY) {
				residGblFitY->setTitle(histTitle);
				_mapFromSensorIDToEfficiencyY.insert(std::make_pair(_sensorIDs.at(i), residGblFitY));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (effGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
void	EUTelProcessorTrackAnalysis::initialiseResidualVsPositionHistograms(){
	int NBinX;
	double MinX;
	double MaxX;
	int NBinY;
	double MinY;
	double MaxY;
	double MinZ;
	double MaxZ;

	std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
    bool isHistoManagerAvailable;
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
	isHistoManagerAvailable = false;


	std::stringstream sstm;
	std::string residGblFitHistName;
	std::string histTitle;
	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
		/////////////////////////////////////////////////////////////////////////////XY residual plots with position
		sstm << "ResidualX" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "ResidualsX. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 200;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 200;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 10;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitX =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitX) {
				residGblFitX->setTitle(histTitle);
				_mapFromSensorIDToHistogramX.insert(std::make_pair(_sensorIDs.at(i), residGblFitX));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "ResidualY" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "ResidualsY. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 200;//every 500 micron there is a bin
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 200;//every 500 micron there is a bin
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 10;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitY =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitY) {
				residGblFitY->setTitle(histTitle);
				_mapFromSensorIDToHistogramY.insert(std::make_pair(_sensorIDs.at(i), residGblFitY));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////The incidence angles for each plane
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceXZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence Global Tx (XZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 180;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IHistogram1D * incidenceGblFitXZ = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX); 

		if (incidenceGblFitXZ){
				incidenceGblFitXZ->setTitle(histTitle);
				_mapFromSensorIDToGloIncXZ.insert(std::make_pair(_sensorIDs.at(i), incidenceGblFitXZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceYZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence Global Ty (YZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 180;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IHistogram1D * incidenceGblFitYZ = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX); 

		if (incidenceGblFitYZ) {
				incidenceGblFitYZ->setTitle(histTitle);
				_mapFromSensorIDToGloIncYZ.insert(std::make_pair(_sensorIDs.at(i), incidenceGblFitYZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	/////////////////////////////////////////////////////////////////////////////////////// 
	/////////////////////////////////////////////////////////////////////////////////////The profile of the p-values with incidence.
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << " Profile of p-values Vs IncidenceXZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Profile of p-values Vs Incidence Angle, global Tx (XZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IProfile1D *pValueVsIncidenceXZ = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D(residGblFitHistName, NBinX, MinX, MaxX, 0, 1); 

		if (pValueVsIncidenceXZ){
				pValueVsIncidenceXZ->setTitle(histTitle);
				_mapFromSensorIDToPValuesVsIncidenceXZ.insert(std::make_pair(_sensorIDs.at(i), pValueVsIncidenceXZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << " Profile of p-values Vs IncidenceYZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Profile of p-values Vs Incidence Angle, global Ty (YZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IProfile1D * pValueVsIncidenceYZ = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D(residGblFitHistName, NBinX, MinX, MaxX, 0,1); 

		if (pValueVsIncidenceYZ) {
				pValueVsIncidenceYZ->setTitle(histTitle);
			_mapFromSensorIDToPValuesVsIncidenceYZ.insert(std::make_pair(_sensorIDs.at(i), pValueVsIncidenceYZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	/////////////////////////////////////////////////////////////////////////////////////// 

	/////////////////////////////////////////////////////////////////////////////////////////p-value with position
	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "P-value vs Position" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "P-value. Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 20;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10.6;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10.6;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 10;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -5.3;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 5.3;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : 0;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 1;
		AIDA::IProfile2D *  pValueHisto =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (pValueHisto) {
				pValueHisto->setTitle(histTitle);
				_mapFromSensorIDToPValueHisto.insert(std::make_pair(_sensorIDs.at(i), pValueHisto));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////Beam Energy
	_beamEnergy  = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("BeamEnergy", 1000, 0, 6); 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////P-value with energy
	_pValueVsBeamEnergy = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D("p-Value Vs Beam Energy", 50, 0, 6, 0,1); 

}
void	EUTelProcessorTrackAnalysis::initialiseGeoInfo(){
  for (size_t i = 0; i < _sensorIDs.size() ; ++i){
    streamlog_out (MESSAGE5)<<"sensor id = "<<i<<" _sensorIDs.at(i) = "<<_sensorIDs.at(i)<<std::endl;
    //get the geoemtry description of the plane
    geo::EUTelGenericPixGeoDescr* geoDescr = geo::gGeometry().getPixGeoDescr( _sensorIDs.at(i) ) ;
    //get the max/min pixel indices
    int minX, minY, maxX, maxY;
    minX = minY = maxX = maxY = 0;
    geoDescr->getPixelIndexRange( minX, maxX, minY, maxY );
    float sizeX, sizeY, sizeZ;
    sizeX= 0.; sizeY=0.; sizeZ=0.;
    geoDescr->getSensitiveSize(sizeX, sizeY, sizeZ);
    streamlog_out (MESSAGE5)<<"sensor id = "<<i<<" has minX, maxX, minY, maxY = "<<minX<<","<< maxX<<","<<  minY<<", "<< maxY <<" and sensitive size = "<<sizeX<<","<< sizeY<<","<<  sizeZ<<std::endl;
    float pitchX =0.0;
    float pitchY =0.0;
    if(_sensorIDs.at(i) >=20){
      sizeX=int(sizeX);

    }//hack by helen.   why is size =20.3????
    streamlog_out (MESSAGE5)<<"sensor id("<<i<<") = "<<_sensorIDs.at(i)<<" has minX, maxX, minY, maxY = "<<minX<<","<< maxX<<","<<  minY<<", "<< maxY <<" and sensitive size = "<<sizeX<<","<< sizeY<<","<<  sizeZ<<std::endl;
    pitchX = sizeX/(maxX-minX+1);
    pitchY = sizeY/(maxY-minY+1);
    streamlog_out (MESSAGE5)<<"pitchX = "<<pitchX<<"   pitchY = "<<pitchY<<std::endl; 
    _mapFromSensorPitchX.insert(std::make_pair(_sensorIDs.at(i), pitchX));
    _mapFromSensorPitchY.insert(std::make_pair(_sensorIDs.at(i), pitchY));
  }
}


void	EUTelProcessorTrackAnalysis::initialiseHitMapHistograms(){
	int NBinX;
	double MinX;
	double MaxX;
	int NBinY;
	double MinY;
	double MaxY;
	double MinZ;
	double MaxZ;

	std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
    bool isHistoManagerAvailable;
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
	isHistoManagerAvailable = false;


	std::stringstream sstm;
	std::string residGblFitHistName;
	std::string histTitle;
	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
		/////////////////////////////////////////////////////////////////////////////XY residual plots with position
		sstm << "HitMap" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "HitMap. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 500;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 500;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -10;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 10;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IHistogram2D *  HitMap =	marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY);
		if (HitMap) {
				HitMap->setTitle(histTitle);
				_mapFromSensorHitMap.insert(std::make_pair(_sensorIDs.at(i), HitMap));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}


}
