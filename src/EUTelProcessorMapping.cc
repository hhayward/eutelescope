/*
 *   Most of this code is based on the EUTelGeometricClusteringProcessor by
 *   Tobias Bisanz
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//eutelescope includes
#include "EUTelProcessorMapping.h"

#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelHistogramManager.h"

//eutel data specific
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"

//eutel geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

//ROOT includes
#include "TGeoShape.h"
#include "TGeoBBox.h"

//marlin includes
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"

//lcio includes
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
//aida includes
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

//system includes
#include <string>
#include <vector>
#include <memory>
//#include <iostream>
#include <cmath>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorMapping::EUTelProcessorMapping(): 
  Processor("EUTelProcessorMapping"), 
  _zsDataCollectionName(""),
  _unmappedCollectionName(""),
  _maskedCollectionName(""),
  _mappedCollectionName(""),
  _maskedOutCollectionName(""),
  _initialUnmappedCollectionSize(0),
  _initialMaskedCollectionSize(0),
  _initialMappedCollectionSize(0),
  _initialMaskedOutCollectionSize(0),
  _iRun(0),
  _iEvt(0),
  _maskChannels(false),
  _fillHistos(false),
  _doLargeClusterMaps(false),
  _doMaskedOutPlots(false),
  _histoInfoFileName(""),
  _noOfDetector(0),
  _ExcludedPlanes(),
  _unmappedHistos(),
  _unmappedSelectedEventsHistos(),
  _maskedOutSpectrumHistos(),
  _maskedHistos(),
  _mappedHistos(),
  _maskedOutHistos(),
  _mappedSelectedEventsHistos(),
  _maskedOutSelectedEventsHistos(),
  _eventMultiplicityHistos(),
  _sensorIDVec(),
  _zsInputDataCollectionVec(NULL),
  _unmappedCollectionVec(NULL),
  _maskedCollectionVec(NULL),
  _mappedCollectionVec(NULL),
  _maskedOutCollectionVec(NULL)
 {
  
  // modify processor description
  _description = "EUTelProcessorGeometricClustering is looking for clusters into a calibrated pixel matrix.";

  // register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName", "Input of Zero Suppressed data",
                           _zsDataCollectionName, std::string ("zsdata") );
     
     // register the output collections: one mapped and one as original
  registerOutputCollection(LCIO::TRACKERDATA, "UnmappedCollectionName", "Unmapped (output) collection name",
                           _unmappedCollectionName, std::string("UnmappedHits"));

  registerOutputCollection(LCIO::TRACKERDATA, "MaskedCollectionName", "Masked (output) collection name",
                           _maskedCollectionName, std::string("maskedHits"));

  registerOutputCollection(LCIO::TRACKERDATA, "MappedCollectionName", "Mapped (output) collection name",
                           _mappedCollectionName, std::string("mappedHits"));

  registerOutputCollection(LCIO::TRACKERDATA, "MaskedOutCollectionName", "MaskedOut (output) collection name",
                           _maskedOutCollectionName, std::string("maskedOutHits"));

  // now the optional parameters
  registerOptionalParameter("DUTid", "The list of sensor ids for mapping",
                             _DUTid, std::vector<int> () );
 // nPixX & nPixY will be used to identify mapping algorithm
  registerOptionalParameter("nPixX", "The list of sensor X dimensions (pixels).",
                             _nPixX, std::vector<int> () );

  registerOptionalParameter("nPixY", "The list of sensor Y dimensions (pixels).",
                             _nPixY, std::vector<int> () );

  registerProcessorParameter("TCut","Time cut in time units of your sensor",
                             _cutT, static_cast<float > ( std::numeric_limits<float>::max() ));

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, std::string( "histoinfo.xml" ) );

  registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",
    _fillHistos, static_cast< bool > ( true ) );

  registerProcessorParameter("MaskingPlots","Switch on or off the channel masking plots",
    _maskChannels, static_cast< bool > ( false ) );

  registerOptionalParameter("ExcludedPlanes", "The list of sensor ids that have to be excluded from the clustering.",
                             _ExcludedPlanes, std::vector<int> () );

  registerOptionalParameter("LargeClusterPlots", "Switch on or off maps of large (>100 pixels) cluster events.",
                             _doLargeClusterMaps, static_cast< bool > ( false ) );

  registerOptionalParameter("UnusedChannelPlots", "Switch on or off plots for unused pixels (from 167x125um**2)",
                             _doMaskedOutPlots, static_cast< bool > ( false ) );

  		_isFirstEvent = true;
}

void EUTelProcessorMapping::init() {
	//this method is called only once even when the rewind is active, it is usually a good idea to
	printParameters ();

	//init new geometry
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,true);

	//set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;

	//the geometry is not yet initialized, so set the corresponding switch to false
	_isGeometryReady = false;
}

void EUTelProcessorMapping::processRunHeader (LCRunHeader * rdr) {

	std::auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );
	runHeader->addProcessor( type() );
  	//increment the run counter
	++_iRun;
}

void EUTelProcessorMapping::initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException ) {

	//set the total number of detector to zero. This number can be different from the one written in the gear description because
	//the input collection can contain only a fraction of all the sensors.

	_noOfDetector = 0;
	_sensorIDVec.clear();

	streamlog_out( DEBUG5 ) << "Initializing geometry" << std::endl;

	for(int d=0; d<_DUTid.size(); d++){
	  _totPixelMap.insert( std::make_pair( _DUTid.at(d), 0) );
	  _totPixelCheck.insert( std::make_pair( _DUTid.at(d), 0) );
	  _foundSelectedEvents.insert( std::make_pair( _DUTid.at(d), 0) );
	  _foundMaskedOutPixels.insert( std::make_pair( _DUTid.at(d), 0) );
	  _noOfDetector++;
	}

  	try 
  	{
		_zsInputDataCollectionVec = dynamic_cast<LCCollectionVec*>( event->getCollection(_zsDataCollectionName) );
		//_noOfDetector += _zsInputDataCollectionVec->getNumberOfElements();
		CellIDDecoder<TrackerDataImpl> cellDecoder(_zsInputDataCollectionVec);

		for ( size_t i = 0; i < _zsInputDataCollectionVec->size(); ++i ) 
		{
			TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( _zsInputDataCollectionVec->getElementAt( i ) ) ;
			_sensorIDVec.push_back( cellDecoder(data)["sensorID"] );
		}
	} 

	catch ( lcio::DataNotAvailableException ) 
	{
		streamlog_out( DEBUG5 ) << "Could not find the input collection: " << _zsDataCollectionName.c_str() << " !" << std::endl;
		return;
	}

    _isGeometryReady = true;
}

void EUTelProcessorMapping::modifyEvent( LCEvent * /* event */ )
{
	return;
}

void EUTelProcessorMapping::readCollections(LCEvent* event)
{
	try 
	{
		_zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
        	streamlog_out ( DEBUG4 ) << "_zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << std::endl;
	} 
	catch ( lcio::DataNotAvailableException )   // do nothing
	{
		streamlog_out ( DEBUG4 ) << "_zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << std::endl;
	}

	try 
	{
		event->getCollection( _zsDataCollectionName ) ;
	} 
	catch (lcio::DataNotAvailableException& e ) 
	{
		streamlog_out(MESSAGE2) << "The current event doesn't contain nZS data collections: skip # " << event->getEventNumber() << std::endl;
		throw SkipEventException(this);
	}
}

void EUTelProcessorMapping::processEvent (LCEvent * event) 
{
	//increment event counter
	++_iEvt;

	//first of all we need to be sure that the geometry is properly initialized!
	if ( !_isGeometryReady ) 
	{
		initializeGeometry( event );
	}

	//read zsInputData collection
	readCollections(event);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	// book the histograms now
	if ( _fillHistos && isFirstEvent() ) 
	{
		bookHistos();
	}
#endif

	EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);
	if ( evt->getEventType() == kEORE ) 
  	{
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  std::endl;
		return;
	}
	else if ( evt->getEventType() == kUNKNOWN ) 
	{
		streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	// prepare a collection of mapped hits can be either a new collection or already existing in the event
	LCCollectionVec* unmappedCollection;
	LCCollectionVec* maskedCollection;
	LCCollectionVec* mappedCollection;
	LCCollectionVec* maskedOutCollection;
	bool unmappedCollectionExists = false;
	bool maskedCollectionExists = false;
	bool mappedCollectionExists = false;
	bool maskedOutCollectionExists = false;
	_initialUnmappedCollectionSize = 0;
	_initialMaskedCollectionSize = 0;
	_initialMappedCollectionSize = 0;
	_initialMaskedOutCollectionSize = 0;
	try 
	{
		unmappedCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _unmappedCollectionName ) );
		unmappedCollectionExists = true;
		_initialUnmappedCollectionSize = unmappedCollection->size();
	} 
	catch ( lcio::DataNotAvailableException& e ) 
	{
		unmappedCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	}

	try 
	{
		maskedCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _maskedCollectionName ) );
		maskedCollectionExists = true;
		_initialMaskedCollectionSize = maskedCollection->size();
	} 
	catch ( lcio::DataNotAvailableException& e ) 
	{
	  maskedCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	}

	try 
	{
		mappedCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _mappedCollectionName ) );
		mappedCollectionExists = true;
		_initialMappedCollectionSize = mappedCollection->size();
	} 
	catch ( lcio::DataNotAvailableException& e ) 
	{
	  mappedCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	}

	try 
	{
		maskedOutCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _maskedOutCollectionName ) );
		maskedOutCollectionExists = true;
		_initialMappedCollectionSize = maskedOutCollection->size();
	} 
	catch ( lcio::DataNotAvailableException& e ) 
	{
	  maskedOutCollection = new LCCollectionVec(LCIO::TRACKERDATA);
	}
  

	//HERE WE ACTUALLY CALL THE MASKING ROUTINE:
	if ( _maskChannels ) {
	  for(int d=0; d<_DUTid.size();d++){
	    hitMasking(evt, maskedCollection, _DUTid.at(d));
	  }
	}

	//HERE WE ACTUALLY CALL THE MAPPING ROUTINE:
	for(int d=0; d<_DUTid.size();d++){
	  hitMapping(evt, unmappedCollection, mappedCollection, maskedOutCollection, _DUTid.at(d));
	}

	// if the dataCollection is not empty add it to the event
	if ( ! unmappedCollectionExists && ( unmappedCollection->size() != _initialUnmappedCollectionSize )) 
	{		evt->addCollection( unmappedCollection, _unmappedCollectionName );	}
	if ( ! maskedCollectionExists && ( maskedCollection->size() != _initialMaskedCollectionSize )) 
	{		evt->addCollection( maskedCollection, _maskedCollectionName );	}
	if ( ! mappedCollectionExists && ( mappedCollection->size() != _initialMappedCollectionSize )) 
	{		evt->addCollection( mappedCollection, _mappedCollectionName );	}
	if ( ! maskedOutCollectionExists && ( maskedOutCollection->size() != _initialMaskedOutCollectionSize )) 
	{		evt->addCollection( maskedOutCollection, _maskedOutCollectionName );	}
  
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	fillHistos(evt);
#endif

	if ( ! unmappedCollectionExists && ( unmappedCollection->size() == _initialUnmappedCollectionSize ) ) 
	{		delete unmappedCollection;	}
	if ( ! maskedCollectionExists && ( maskedCollection->size() == _initialMaskedCollectionSize ) ) 
	{		delete maskedCollection;	}
	if ( ! mappedCollectionExists && ( mappedCollection->size() == _initialMappedCollectionSize ) ) 
	{		delete mappedCollection;	}
	if ( ! maskedOutCollectionExists && ( maskedOutCollection->size() == _initialMaskedOutCollectionSize ) ) 
	{		delete maskedOutCollection;	}

	_isFirstEvent = false;
}

void EUTelProcessorMapping::hitMasking( LCEvent * evt , LCCollectionVec * maskedCollection, int DUTid) 
{
	// prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder( _zsInputDataCollectionVec );

  bool isDummyAlreadyExisting = false;
  LCCollectionVec* sparsePixelCollectionVec = NULL;
  
  try 
    {
      sparsePixelCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
      isDummyAlreadyExisting = true ;
    } 
  catch (lcio::DataNotAvailableException& e) 
    {
      sparsePixelCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
      isDummyAlreadyExisting = false;
    }
  
  for ( unsigned int idetector = 0 ; idetector < _zsInputDataCollectionVec->size(); idetector++ ) 
    {
      // get the TrackerData and guess which kind of sparsified data it contains.
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( _zsInputDataCollectionVec->getElementAt( idetector ) );
      SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
      int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

      //if this is an excluded sensor go to the next element
      bool foundexcludedsensor = false;
      for(size_t iexclude = 0; iexclude < _ExcludedPlanes.size(); ++iexclude)
	{
	  if(_ExcludedPlanes[iexclude] == sensorID)
	    {
	      foundexcludedsensor = true;
	    }
	}
      
      if(foundexcludedsensor)       
	{		
	  continue;
	}
      
      int nRow=-1, nCol=-1;
      bool foundDUT=false;
      //find mapping by identifying DUT number and hence corrsponding Pix array (by ordered input vectors)
      for(size_t id= 0; id < _DUTid.size(); id++){
	if(_DUTid[id]==DUTid){
	  foundDUT=true;
	  nRow=_nPixY[id];
	  nCol=_nPixX[id];
	  streamlog_out ( DEBUG2 ) << "identified DUTid "<<DUTid<<" with  _DUTid " << std::endl;
	}
      }
      if(!foundDUT){ nRow=336; nCol=80; 
	streamlog_out ( DEBUG2 ) << "could not identify DUTid " << DUTid << std::endl; 
      }
      //check with sensorID that this is appropriate
      if(sensorID%10!=DUTid%10){
	  foundDUT=false;
	  streamlog_out ( DEBUG2 ) << "Sensor id " << sensorID << " does not match DUTid "<<DUTid<< std::endl; 
      }
      if(!foundDUT){
	streamlog_out ( DEBUG2 ) << "Skipping sensor id " << sensorID << std::endl; 
	continue;
      }
      
      if ( type == kEUTelGenericSparsePixel ) 
	{
	  
	  // now prepare the EUTelescope interface to sparsified data.
	  std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );
	  
	  streamlog_out ( DEBUG2 ) << "Processing sparse data on detector " << sensorID << " with " << sparseData->size() << " pixels " << std::endl;
	  
	  int hitPixelsInEvent = sparseData->size();
	  
	  //	create cell encoders to set sensorID and pixel type
	  CellIDEncoder< TrackerDataImpl > zsMaskedDataEncoder   ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, maskedCollection  );
	  // prepare a new TrackerData object for the ZS data
	  // it contains all the hits for a particular sensor in one event
	  std::auto_ptr<lcio::TrackerDataImpl > zsMaskedFrame( new lcio::TrackerDataImpl );
	  zsMaskedDataEncoder["sensorID"] = DUTid;
	  zsMaskedDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
	  // set some values of "Cells" for this object
	  zsMaskedDataEncoder.setCellID( zsMaskedFrame.get() );
	  std::auto_ptr< eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > >
	    sparseMaskedFrame( new eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > ( zsMaskedFrame.get() ) );
	  EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;
	  
	  //This for-loop loads all the hits of the given event and detector plane and stores them 
	  for(int i = 0; i < hitPixelsInEvent; ++i )
	    {
	      //Load the information of the hit pixel into genericPixel
	      sparseData->getSparsePixelAt( i, genericPixel );		
	      
	      //Now we need to process the pixel
	      if ( genericPixel ) 
		{

		  // get the cluster size in X and Y separately and plot it:
		  int xPos= genericPixel->getXCoord(); int yPos=genericPixel->getYCoord();
		  int row = yPos, col = xPos;
		  // MASKING SECTION//



		  if(nCol==120 && nRow==136){ //example mapping for 167x125 device (not universally applicable!)
		    
		    if(DUTid==30){ // specific to CPV 167x125 device readout (offset 1)
		      if( (col%4==1 || col%4==2) && (row%5==2 || row%5==0) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==3 || row%5==0) ){ continue; } //no 2s or 4s
		    }
		    else if(DUTid==40){ // specific to CPV 167x125 device readout (offset 2)
		      if( (col%4==1 || col%4==2) && (row%5==3 || row%5==1) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==4 || row%5==1) ){ continue; } //no 2s or 4s
		    }
		    else if(DUTid==50){ // specific to CPV 167x125 device readout (offset 3)
		      if( (col%4==1 || col%4==2) && (row%5==4 || row%5==2) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==0 || row%5==2) ){ continue; } //no 2s or 4s
		    }	  
		    else if(DUTid==60||DUTid==62){ // specific to CPV 167x125 device readout (offset 4)
		      if( (col%4==1 || col%4==2) && (row%5==0 || row%5==3) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==1 || row%5==3) ){ continue; } //no 2s or 4s
		    }	  
		    else if(DUTid==70||DUTid==72){ // specific to CPV 167x125 device readout (offset 4)
		      if( (col%4==1 || col%4==2) && (row%5==0 || row%5==3) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==1 || row%5==3) ){ continue; } //no 2s or 4s
		    }	  
		    else if(DUTid==20||DUTid==22){ // specific to CPV 167x125 device readout (offset 4)
		      if( (col%4==1 || col%4==2) && (row%5==0 || row%5==3) ){ continue; } //no 6s or 9s
		      if( (col%4==0 || col%4==3) && (row%5==1 || row%5==3) ){ continue; } //no 2s or 4s
		    }
		  }
		  // add editted pixel to mapped pixel map
		  sparseMaskedFrame->addSparsePixel( genericPixel );
		  
		  streamlog_out ( DEBUG2 ) << "storing masked pixels ("<<DUTid <<")" << genericPixel->getXCoord() <<", " << genericPixel->getYCoord() << std::endl;

		}
	    }
	  
	  // write TrackerData object that contains info from one sensor to LCIO collection
	  maskedCollection->push_back( zsMaskedFrame.release() );
	  
	}
    } 
} 

void EUTelProcessorMapping::hitMapping( LCEvent * evt , LCCollectionVec * unmappedCollection, LCCollectionVec * mappedCollection, LCCollectionVec * maskedOutCollection, int DUTid) 
{
	// prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder( _zsInputDataCollectionVec );

  bool isDummyAlreadyExisting = false;
  LCCollectionVec* sparsePixelCollectionVec = NULL;
  
  try 
    {
      sparsePixelCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
      isDummyAlreadyExisting = true ;
    } 
  catch (lcio::DataNotAvailableException& e) 
    {
      sparsePixelCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
      isDummyAlreadyExisting = false;
    }
  
  for ( unsigned int idetector = 0 ; idetector < _zsInputDataCollectionVec->size(); idetector++ ) 
    {
      // get the TrackerData and guess which kind of sparsified data it contains.
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( _zsInputDataCollectionVec->getElementAt( idetector ) );
      SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
      int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

      //if this is an excluded sensor go to the next element
      bool foundexcludedsensor = false;
      for(size_t iexclude = 0; iexclude < _ExcludedPlanes.size(); ++iexclude)
	{
	  if(_ExcludedPlanes[iexclude] == sensorID)
	    {
	      foundexcludedsensor = true;
	    }
	}
      
      if(foundexcludedsensor)       
	{		
	  continue;
	}
      
      int nRow=-1, nCol=-1;
      bool foundDUT=false;
      //find mapping by identifying DUT number and hence corrsponding Pix array (by ordered input vectors)
      for(size_t id= 0; id < _DUTid.size(); id++){
	if(_DUTid[id]==DUTid){
	  foundDUT=true;
	  nRow=_nPixY[id];
	  nCol=_nPixX[id];
	  streamlog_out ( DEBUG2 ) << "identified DUTid "<<DUTid<<" with  _DUTid " << _DUTid[id] << ": " <<nCol<<"x"<<nRow<< std::endl;
	}
      }
      if(!foundDUT){ nRow=336; nCol=80; 
	streamlog_out ( DEBUG2 ) << "could not identify DUTid " << DUTid << ". Standard geometry used"<< std::endl; 
      }
      //check with sensorID that this is appropriate
      if(sensorID%10!=DUTid%10){
	  foundDUT=false;
	  streamlog_out ( DEBUG2 ) << "Sensor id " << sensorID << " does not match DUTid "<<DUTid<< std::endl; 
      }
      if(!foundDUT){
	streamlog_out ( DEBUG2 ) << "Skipping sensor id " << sensorID << std::endl; 
	continue;
      }
      
      if ( type == kEUTelGenericSparsePixel ) 
	{
	  
	  // now prepare the EUTelescope interface to sparsified data.
	  std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );
	  
	  streamlog_out ( DEBUG2 ) << "Processing sparse data on detector " << sensorID << " with " << sparseData->size() << " pixels " << std::endl;
	  
	  int hitPixelsInEvent = sparseData->size();
	  
	  //	create cell encoders to set sensorID and pixel type
	  CellIDEncoder< TrackerDataImpl > zsUnmappedDataEncoder   ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, unmappedCollection  );
	  std::auto_ptr<lcio::TrackerDataImpl > zsUnmappedFrame( new lcio::TrackerDataImpl );
	  zsUnmappedDataEncoder["sensorID"] = DUTid;
	  zsUnmappedDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
	  zsUnmappedDataEncoder.setCellID( zsUnmappedFrame.get() );
	  //	create cell encoders to set sensorID and pixel type
	  std::auto_ptr< eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > >
	    sparseUnmappedFrame( new eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > ( zsUnmappedFrame.get() ) );

	  CellIDEncoder< TrackerDataImpl > zsMappedDataEncoder   ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, mappedCollection  );
	  // prepare a new TrackerData object for the ZS data
	  // it contains all the hits for a particular sensor in one event
	  std::auto_ptr<lcio::TrackerDataImpl > zsMappedFrame( new lcio::TrackerDataImpl );
	  zsMappedDataEncoder["sensorID"] = DUTid;
	  zsMappedDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
	  // set some values of "Cells" for this object
	  zsMappedDataEncoder.setCellID( zsMappedFrame.get() );
	  std::auto_ptr< eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > >
	    sparseMappedFrame( new eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > ( zsMappedFrame.get() ) );

	  CellIDEncoder< TrackerDataImpl > zsMaskedOutDataEncoder   ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, maskedOutCollection  );
	  // prepare a new TrackerData object for the ZS data
	  // it contains all the hits for a particular sensor in one event
	  std::auto_ptr<lcio::TrackerDataImpl > zsMaskedOutFrame( new lcio::TrackerDataImpl );
	  zsMaskedOutDataEncoder["sensorID"] = DUTid;
	  zsMaskedOutDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
	  // set some values of "Cells" for this object
	  zsMaskedOutDataEncoder.setCellID( zsMaskedOutFrame.get() );
	  std::auto_ptr< eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > >
	    sparseMaskedOutFrame( new eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > ( zsMaskedOutFrame.get() ) );

	  EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;
	  EUTelGenericSparsePixel* genericMappedPixel = new EUTelGenericSparsePixel;
	  
	  //This for-loop loads all the hits of the given event and detector plane and stores them 
	  for(int i = 0; i < hitPixelsInEvent; ++i )
	    {
	      //Load the information of the hit pixel into genericPixel
	      sparseData->getSparsePixelAt( i, genericPixel );		
	      
	      //Now we need to process the pixel
	      if ( genericPixel ) 
		{
		  //save the data to unmapped pixel map
		  sparseUnmappedFrame->addSparsePixel( genericPixel );

		  // get the cluster size in X and Y separately and plot it:
		  int xPos= genericPixel->getXCoord(); int yPos=genericPixel->getYCoord();
		  int row = yPos, col = xPos;
		  // MAPPING SECTION//
		  int geoRow=-1, geoCol=-1;
		  bool maskedOut=false;
		  
		  if(nCol==40 && nRow==672){ //mapping for UK 500x25 device (not universally applicable!)
		    if(col%2==0){
		      geoRow=row*2/1 +1;
		      geoCol=(col)*1/2;
		    }
		    else{
		      geoRow=row*2/1;
		      geoCol=(col-1)*1/2; 
		    }
		  }
		  else if(nCol==120 && (nRow==136 ||nRow==134 )){ //mapping for UK 167x125 device (not universally applicable!)
		    // specific to CPV 167x125 device readout (offset 4)
		    if( (col%4==1 || col%4==2) && (row%5==0 || row%5==3) ){ maskedOut=true; } //no 6s or 9s
		    if( (col%4==0 || col%4==3) && (row%5==1 || row%5==3) ){ maskedOut=true; } //no 2s or 4s
		    if(!maskedOut){
		      if(col%4==0 || col%4==3){ // 0-4 columns
			if(row%5==4){ //row 0s
			  geoRow=(row*2 +2)/5;
			  if(col%4==0){
			    geoCol=col*3/2;
			    }
			  else if(col%4==3){
			    geoCol=(col*3 +1)/2;
			  }
			}
			else if(row%5==0){ //row 1s (drop 2)
			  geoRow=(row*2)/5 ;
			  if(col%4==0){
			    geoCol=(col*3 +2)/2;
			  }
			  else if(col%4==3){
			    geoCol=(col*3 -1)/2;
			  }
			}
			else if(row%5==2){ // row 3s
			  geoRow=(row*2 +1)/5 ;
			  if(col%4==0){
			    geoCol=col*3/2;
			  }
			  else if(col%4==3){
			    geoCol=(col*3 +1)/2;
			    }
			}
		      }
		      if(col%4==1 || col%4==2){ // 5-9 columns
			if(row%5==4){ //row 5s (drop 1)
			  geoRow=(row*2 +2)/5 ;
			  if(col%4==1){
			    geoCol=(col*3 +1)/2;
			  }
			  else if(col%4==2){
			    geoCol=col*3/2;
			  }
			}
			else if(row%5==1){ //row 7s (drop 2)
			  geoRow=(row*2 +3)/5 ;
			  if(col%4==1){
			    geoCol=(col*3 -1)/2;
			  }
			  else if(col%4==2){
			    geoCol=(col*3 +2)/2;
			  }
			}
			else if(row%5==2){ //row 8s (drop 3)
			  geoRow=(row*2 +1)/5 ;
			  if(col%4==1){
			    geoCol=(col*3 +1)/2;
			  }
			  else if(col%4==2){
			    geoCol=col*3/2;
			  }
			}
		      }
		    }
		  }	
		  else if(nCol==160 && nRow==168){ //mapping for UK 125x100 device (not universally applicable!)
		    if(row%4==0){
		      geoRow=row*1/2;
		      if(col%2==0){ geoCol=col*2/1 ; }
		      if(col%2==1){ geoCol=col*2/1 +1; }
		    }
		    if(row%4==1){
		      geoRow=(row-1)*1/2;
		      if(col%2==0){ geoCol=col*2/1 +1; }
		      if(col%2==1){ geoCol=col*2/1 ; }
		    }
		    if(row%4==2){
		      geoRow=(row)*1/2;
		      if(col%2==0){ geoCol=col*2/1 +1; }
		      if(col%2==1){ geoCol=col*2/1 ; }
		    }
		    if(row%4==3){
		      geoRow=(row-1)*1/2;
		      if(col%2==0){ geoCol=col*2/1 ; }
		      if(col%2==1){ geoCol=col*2/1 +1; }
		    }
		  }
		  else{ geoRow=row; geoCol=col; }
		  
		  if(maskedOut ){ geoRow=row; geoCol=col; }
		  // get data pixel for mapping coordinates
		    sparseData->getSparsePixelAt( i, genericMappedPixel );	
		    genericMappedPixel->setXCoord(geoCol);
		    genericMappedPixel->setYCoord(geoRow);
		    // add editted pixel to mapped pixel map
		    
		    if(maskedOut ){
		      sparseMaskedOutFrame->addSparsePixel( genericMappedPixel );	  
		      streamlog_out ( DEBUG2 ) << "storing maskedout pixels ("<<DUTid <<")" << genericMappedPixel->getXCoord() <<", " << genericMappedPixel->getYCoord() << std::endl;
		_foundMaskedOutPixels[DUTid] += 1;
		    }
		    else{
		      sparseMappedFrame->addSparsePixel( genericMappedPixel );	  
		      streamlog_out ( DEBUG2 ) << "storing mapped pixels ("<<DUTid <<")" << genericMappedPixel->getXCoord() <<", " << genericMappedPixel->getYCoord() << std::endl;
		    }
		    
		    streamlog_out ( DEBUG2 ) << "storing unmapped pixels ("<<DUTid <<") " << genericPixel->getXCoord() <<", " << genericPixel->getYCoord() << std::endl;
		    
		    // last but not least increment the totClusterMap
		    _totPixelMap[ DUTid ] += 1;
		}
	    }
	  
	  // write TrackerData object that contains info from one sensor to LCIO collection
	  unmappedCollection->push_back( zsUnmappedFrame.release() );
	  mappedCollection->push_back( zsMappedFrame.release() );
	  maskedOutCollection->push_back( zsMaskedOutFrame.release() );
	  
	}
    } 
} 

void EUTelProcessorMapping::check (LCEvent * /* evt */) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}

void EUTelProcessorMapping::end() {
  
	std::map<int,int>::iterator iter = _totPixelMap.begin();
	while( iter != _totPixelMap.end() )
	{
		streamlog_out ( MESSAGE4 ) << "Found " << iter->second << " pixels on detector " << iter->first << std::endl;
		++iter;
	}
	iter = _totPixelCheck.begin();
	while( iter != _totPixelCheck.end() )
	{
		streamlog_out ( MESSAGE4 ) << "Found " << iter->second << " check on detector " << iter->first << std::endl;
		++iter;
	}
	iter = _foundSelectedEvents.begin();
	while( iter != _foundSelectedEvents.end() )
	{
		streamlog_out ( MESSAGE4 ) << "Found " << iter->second << " selected events on detector " << iter->first << std::endl;
		++iter;
	}
	if(_doMaskedOutPlots){
	  iter = _foundMaskedOutPixels.begin();
	  while( iter != _foundMaskedOutPixels.end() )
	    {
	      streamlog_out ( MESSAGE4 ) << "Found " << iter->second << " masked out pixels on detector " << iter->first << std::endl;
		++iter;
	    }
	}
	streamlog_out ( MESSAGE4 ) <<  "Successfully finished" << std::endl;
}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelProcessorMapping::fillHistos (LCEvent * evt) 
{ //plot hitmaps to check everything has worked out
  EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type = eutelEvent->getEventType();
  

  if ( type == kEORE ) 
    {
      streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << std::endl;
      return;
    }
  else if ( type == kUNKNOWN ) 
    {
      // if it is unknown we had already issued a warning to the user at
      // the beginning of the processEvent. If we get to here, it means
      // that the assumption that the event was a data event was
      // correct, so no harm to continue...
    }
  
  bool selected=false;// make true if >100 hit pixels in event

  try 
    { // UNMAPPED plots
      LCCollectionVec* _unmappedCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_unmappedCollectionName));
      CellIDDecoder<TrackerDataImpl > cellDecoder(_unmappedCollectionVec);
      
      std::map<int, int> eventCounterMap;
      
	for ( unsigned int idetector = 0 ; idetector < _DUTid.size(); idetector++ ) 
	{
	  TrackerDataImpl* hits = dynamic_cast<TrackerDataImpl*> ( _unmappedCollectionVec->getElementAt(idetector) );
	  SparsePixelType type  = static_cast<SparsePixelType> ( static_cast<int> ( cellDecoder(hits)["sparsePixelType"] ));
	  int detectorID = static_cast<int>( cellDecoder(hits)["sensorID"] );
	  	  
	  if( type == kEUTelGenericSparsePixel ) {

	    //if this key doesn't exist yet it will be value initialized, this is desired, for int this is 0!
	    eventCounterMap[detectorID]++;
	    
	    //if this is an excluded sensor go to the next element
	    bool foundexcludedsensor = false;
	    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
	      {
		if(_ExcludedPlanes[i] == detectorID)
		  {
		    foundexcludedsensor = true;
		    break;
		  }
	      }
	    if(foundexcludedsensor){ continue; }
	    // now prepare the EUTelescope interface to sparsified data.
	    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > hitData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hits ) );
	    streamlog_out ( DEBUG2 ) << "Processing hits for unmapping on detector " << detectorID << " with " << hitData->size() << " pixels " << std::endl;
	    
	    int hitPixelsInEvent = hitData->size();
	    std::vector<EUTelGeometricPixel> hitPixelVec;
	    EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;

	    /*for( int s=0; s<_selectedEvents[detectorID].size(); s++){
	      if(evt->getEventNumber()== _selectedEvents[detectorID].at(s)){ selected=true; 
		_foundSelectedEvents[detectorID] += 1; break; }
	      }*/
	    
	    	//fill the event multiplicity here
	    std::string tempHistoName;
	    AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _eventMultiplicityHistos[detectorID] );
	    if ( histo ) { histo->fill( hitPixelsInEvent ); }

	    if(hitPixelsInEvent>100){ selected=true; 
		_foundSelectedEvents[detectorID] += 1; }
	    
	    //This for-loop loads all the hits of the given event and detector plane and stores them as GeometricPixels
	    for(int i = 0; i < hitPixelsInEvent; ++i )
	      {
		//Load the information of the hit pixel into genericPixel
		hitData->getSparsePixelAt( i, genericPixel );
		EUTelGeometricPixel hitPixel( *genericPixel );
		
		// get the cluster size in X and Y separately and plot it:
		int xPos= genericPixel->getXCoord();
		int yPos= genericPixel->getYCoord();
		streamlog_out ( DEBUG2 ) << "unmapped hits_" <<i<<" : "<< xPos<< ", "<<yPos<<std::endl;
		
		//Do all the plots
		(dynamic_cast<AIDA::IHistogram2D*> (_unmappedHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.);		
		if(_doLargeClusterMaps && selected){
		  (dynamic_cast<AIDA::IHistogram2D*> (_unmappedSelectedEventsHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.); }
	      }
		delete genericPixel;
	  }
	  else 
	    {
	      streamlog_out ( ERROR4 ) <<  "Unknown pixel type:" << type <<std::endl;
	      throw UnknownDataTypeException("Pixel type unknown");
	    }	
	}
    }
  catch (lcio::DataNotAvailableException& e) 
    {
      return;
    }
  
  try 
    { //MAPPED Plots
      LCCollectionVec* _mappedCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_mappedCollectionName));
      CellIDDecoder<TrackerDataImpl > cellDecoder(_mappedCollectionVec);
      
      std::map<int, int> eventCounterMap;
      
	for ( unsigned int idetector = 0 ; idetector < _DUTid.size(); idetector++ ) 
	{
	  TrackerDataImpl* hits = dynamic_cast<TrackerDataImpl*> ( _mappedCollectionVec->getElementAt(idetector) );
	  SparsePixelType type  = static_cast<SparsePixelType> ( static_cast<int> ( cellDecoder(hits)["sparsePixelType"] ));
	  int detectorID = static_cast<int>( cellDecoder(hits)["sensorID"] );
	  	  
	  if( type == kEUTelGenericSparsePixel ) {

	    //if this key doesn't exist yet it will be value initialized, this is desired, for int this is 0!
	    eventCounterMap[detectorID]++;
	    
	    //if this is an excluded sensor go to the next element
	    bool foundexcludedsensor = false;
	    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
	      {
		if(_ExcludedPlanes[i] == detectorID)
		  {
		    foundexcludedsensor = true;
		    break;
		  }
	      }
	    if(foundexcludedsensor){ continue; }
	    // now prepare the EUTelescope interface to sparsified data.
	    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > hitData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hits ) );
	    streamlog_out ( DEBUG2 ) << "Processing hits for mapping on detector " << detectorID << " with " << hitData->size() << " pixels " << std::endl;
	    
	    int hitPixelsInEvent = hitData->size();
	    EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;

	    //bool selected=false;
	    /*for( int s=0; s<_selectedEvents[detectorID].size(); s++){
	      if(evt->getEventNumber()== _selectedEvents[detectorID].at(s)){ selected=true; 
		_foundSelectedEvents[detectorID] += 1; break; }
	      }*/
	    //if(hitPixelsInEvent>100){ selected=true; [detectorID] += 1; }
	    
	    
	    //This for-loop loads all the hits of the given event and detector plane and stores them as GeometricPixels
	    for(int i = 0; i < hitPixelsInEvent; ++i )
	      {
		//Load the information of the hit pixel into genericPixel
		hitData->getSparsePixelAt( i, genericPixel );
		
		// get the cluster size in X and Y separately and plot it:
		int xPos= genericPixel->getXCoord();
		int yPos= genericPixel->getYCoord();
		
		streamlog_out ( DEBUG2 ) << "mapped hits_" <<i<<" : "<< xPos<< ", "<<yPos<<std::endl;
		//Do all the plots
		(dynamic_cast<AIDA::IHistogram2D*> (_mappedHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.);
		if(_doLargeClusterMaps && selected){
		  (dynamic_cast<AIDA::IHistogram2D*> (_mappedSelectedEventsHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.); }
	      }
		delete genericPixel;
	  }
	  else 
	    {
	      streamlog_out ( ERROR4 ) <<  "Unknown pixel type:" << type <<std::endl;
	      throw UnknownDataTypeException("Pixel type unknown");
	    }	
	}

    }
  catch (lcio::DataNotAvailableException& e) 
    {
      return;
    }
  

  if(_doMaskedOutPlots){ 
  try 
    { //MASKEDOUT Plots
      LCCollectionVec* _maskedOutCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_maskedOutCollectionName));
      CellIDDecoder<TrackerDataImpl > cellDecoder(_maskedOutCollectionVec);
      
      std::map<int, int> eventCounterMap;
      
	for ( unsigned int idetector = 0 ; idetector < _DUTid.size(); idetector++ ) 
	{
	  TrackerDataImpl* hits = dynamic_cast<TrackerDataImpl*> ( _maskedOutCollectionVec->getElementAt(idetector) );
	  SparsePixelType type  = static_cast<SparsePixelType> ( static_cast<int> ( cellDecoder(hits)["sparsePixelType"] ));
	  int detectorID = static_cast<int>( cellDecoder(hits)["sensorID"] );
	  	  
	  if( type == kEUTelGenericSparsePixel ) {

	    //if this key doesn't exist yet it will be value initialized, this is desired, for int this is 0!
	    eventCounterMap[detectorID]++;
	    
	    //if this is an excluded sensor go to the next element
	    bool foundexcludedsensor = false;
	    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
	      {
		if(_ExcludedPlanes[i] == detectorID)
		  {
		    foundexcludedsensor = true;
		    break;
		  }
	      }
	    if(foundexcludedsensor){ continue; }
	    // now prepare the EUTelescope interface to sparsified data.
	    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > hitData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hits ) );
	    streamlog_out ( DEBUG2 ) << "Processing hits for mapping on detector " << detectorID << " with " << hitData->size() << " pixels " << std::endl;
	    
	    int hitPixelsInEvent = hitData->size();
	    EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;

	    //bool selected=false;
	    /*for( int s=0; s<_selectedEvents[detectorID].size(); s++){
	      if(evt->getEventNumber()== _selectedEvents[detectorID].at(s)){ selected=true; 
		_foundSelectedEvents[detectorID] += 1; break; }
	      }*/
	    //if(hitPixelsInEvent>100){ selected=true; [detectorID] += 1; }
	    
	    
	    //This for-loop loads all the hits of the given event and detector plane and stores them as GeometricPixels
	    for(int i = 0; i < hitPixelsInEvent; ++i )
	      {
		//Load the information of the hit pixel into genericPixel
		hitData->getSparsePixelAt( i, genericPixel );
		
		// get the cluster size in X and Y separately and plot it:
		int xPos= genericPixel->getXCoord();
		int yPos= genericPixel->getYCoord();
		
		streamlog_out ( DEBUG2 ) << "maskedout hits_" <<i<<" : "<< xPos<< ", "<<yPos<<std::endl;
		//Do all the plots
		_foundMaskedOutPixels[detectorID] += 1;
		(dynamic_cast<AIDA::IHistogram2D*> (_maskedOutHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.);
		(dynamic_cast<AIDA::IHistogram1D*> (_maskedOutSpectrumHistos[detectorID]))->fill(static_cast<double >(genericPixel->getSignal()));
		if(_doLargeClusterMaps && selected){
		  (dynamic_cast<AIDA::IHistogram2D*> (_maskedOutSelectedEventsHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.); }
	      }
		delete genericPixel;
	  }
	  else 
	    {
	      streamlog_out ( ERROR4 ) <<  "Unknown pixel type:" << type <<std::endl;
	      throw UnknownDataTypeException("Pixel type unknown");
	    }	
	}

    }
  catch (lcio::DataNotAvailableException& e) 
    {
      return;
    }
  }

  try 
    { //MASKED Plots
      LCCollectionVec* _maskedCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_maskedCollectionName));
      CellIDDecoder<TrackerDataImpl > cellDecoder(_maskedCollectionVec);
      
      std::map<int, int> eventCounterMap;
      
	for ( unsigned int idetector = 0 ; idetector < _DUTid.size(); idetector++ ) 
	{
	  TrackerDataImpl* hits = dynamic_cast<TrackerDataImpl*> ( _maskedCollectionVec->getElementAt(idetector) );
	  SparsePixelType type  = static_cast<SparsePixelType> ( static_cast<int> ( cellDecoder(hits)["sparsePixelType"] ));
	  int detectorID = static_cast<int>( cellDecoder(hits)["sensorID"] );
	  	  
	  if( type == kEUTelGenericSparsePixel ) {

	    //if this key doesn't exist yet it will be value initialized, this is desired, for int this is 0!
	    eventCounterMap[detectorID]++;
	    
	    //if this is an excluded sensor go to the next element
	    bool foundexcludedsensor = false;
	    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
	      {
		if(_ExcludedPlanes[i] == detectorID)
		  {
		    foundexcludedsensor = true;
		    break;
		  }
	      }
	    if(foundexcludedsensor){ continue; }
	    // now prepare the EUTelescope interface to sparsified data.
	    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > hitData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hits ) );
	    streamlog_out ( DEBUG2 ) << "Processing hits for masking on detector " << detectorID << " with " << hitData->size() << " pixels " << std::endl;
	    
	    int hitPixelsInEvent = hitData->size();
	    EUTelGenericSparsePixel* genericPixel = new EUTelGenericSparsePixel;
	    
	    //This for-loop loads all the hits of the given event and detector plane and stores them as GeometricPixels
	    for(int i = 0; i < hitPixelsInEvent; ++i )
	      {
		//Load the information of the hit pixel into genericPixel
		hitData->getSparsePixelAt( i, genericPixel );
		
		// get the cluster size in X and Y separately and plot it:
		int xPos= genericPixel->getXCoord();
		int yPos= genericPixel->getYCoord();
		
		streamlog_out ( DEBUG2 ) << "mapped hits_" <<i<<" : "<< xPos<< ", "<<yPos<<std::endl;
		//Do all the plots
		(dynamic_cast<AIDA::IHistogram2D*> (_maskedHistos[detectorID]))->fill(static_cast<double >(xPos), static_cast<double >(yPos), 1.);
	      }
		delete genericPixel;
	  }
	  else 
	    {
	      streamlog_out ( ERROR4 ) <<  "Unknown pixel type:" << type <<std::endl;
	      throw UnknownDataTypeException("Pixel type unknown");
	    }	
	}
    }
  catch (lcio::DataNotAvailableException& e) 
    {
      return;
    }
}
#endif

#ifdef MARLIN_USE_AIDA
void EUTelProcessorMapping::bookHistos() {

  // histograms are grouped in loops and detectors
  streamlog_out ( DEBUG5 )  << "Booking histograms " << std::endl;
  std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ) );
  EUTelHistogramInfo* histoInfo;
  bool isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( std::ios::failure& e) {
    streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"
                               << "Continuing without histogram manager"  << std::endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( WARNING2 ) << e.what() << "\n"
                               << "Continuing without histogram manager" << std::endl;
    isHistoManagerAvailable = false;
  }

  // define our histogram names
  std::string _hitNumberHistoName             = "hitNumber";
  std::string _unmappedHistoName             = "ummappedHitMap";
  std::string _unmappedSelectedEventsHistoName             = "ummappedSelectedEventsHitMap";
  std::string _maskedHistoName         = "maskedHitMap";
  std::string _mappedHistoName         = "mappedHitMap";
  std::string _maskedOutHistoName         = "maskedOutHitMap";
  std::string _mappedSelectedEventsHistoName         = "mappedSelectedEventsHitMap";
  std::string _maskedOutSelectedEventsHistoName         = "maskedOutSelectedEventsHitMap";
  std::string _maskedOutSpectrumHistoName         = "maskedOutSpectrum";
  std::string _eventMultiplicityHistoName  = "eventMultiplicity";

  std::string tempHistoName;
  std::string basePath;

	for (size_t iDetector = 0; iDetector < _DUTid.size(); iDetector++) 
	{
		int sensorID = _DUTid.at( iDetector );
        geo::EUTelGenericPixGeoDescr* geoDescr =  ( geo::gGeometry().getPixGeoDescr( sensorID ) );
        
        int nRow=-1, nCol=-1;
        bool foundDUT=false;
        //find mapping
        for(size_t id= 0; id < _DUTid.size(); id++){
            if(sensorID==_DUTid[id]){
                foundDUT=true;
                nRow=_nPixY[id];
                nCol=_nPixX[id];
            }
        }
        if(!foundDUT){ nRow=336; nCol=80; }
        streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... Ncols: "<<nCol<<", Nrows: "<<nRow<< std::endl;

		//now that we know which is the sensorID, we can ask which are the minX, minY, maxX and maxY.
		int minX, minY, maxX, maxY;
		minX = minY = 0; maxX = nCol-1; maxY = nRow-1;
		int uMinX, uMinY, uMaxX, uMaxY;
		uMinX = uMinY = uMaxX = uMaxY = 0;
		 
		geoDescr->getPixelIndexRange( uMinX, uMaxX, uMinY, uMaxY );

		basePath = "detector_" + to_string( sensorID );
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");
        
		streamlog_out ( MESSAGE4 ) <<  "unmapped sensor "<<sensorID <<"... uMaxX: "<<uMaxX<<", uMaxY: "<<uMaxY<< std::endl;
		tempHistoName = _unmappedHistoName + "_d" + to_string( sensorID );
		int     xBin = uMaxX - uMinX + 1;
		double  xMin = static_cast<double >( uMinX ) - 0.5 ;
		double  xMax = static_cast<double >( uMaxX ) + 0.5;
		int     yBin = uMaxY - uMinY + 1;
		double  yMin = static_cast<double >( uMinY ) - 0.5;
        double  yMax = static_cast<double >( uMaxY ) + 0.5;
        //streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... xMax: "<<xMax<<", yMax: "<<yMax<< std::endl;
		AIDA::IHistogram2D * unmappedHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_unmappedHistos.insert(std::make_pair(sensorID, unmappedHisto));
		unmappedHisto->setTitle("UNMAPPED Pixel Index Hit Map;X Index [#];Y Index [#];Count [#]");

		if(_doLargeClusterMaps){
		  tempHistoName = _unmappedSelectedEventsHistoName + "_d" + to_string( sensorID );
		  //streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... xMax: "<<xMax<<", yMax: "<<yMax<< std::endl;
		  AIDA::IHistogram2D * unmappedSelectedEventsHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		  _unmappedSelectedEventsHistos.insert(std::make_pair(sensorID, unmappedSelectedEventsHisto));
		  unmappedSelectedEventsHisto->setTitle("UNMAPPED Pixel Index Hit Map for Selected Events ;X Index [#];Y Index [#];Count [#]");
		}
        
		if(_maskChannels){
		  streamlog_out ( MESSAGE4 ) <<  "masked sensor "<<sensorID <<"... uMaxX: "<<uMaxX<<", uMaxY: "<<uMaxY<< std::endl;
		  tempHistoName = _maskedHistoName + "_d" + to_string( sensorID );
		  AIDA::IHistogram2D * maskedHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		  _maskedHistos.insert(std::make_pair(sensorID, maskedHisto));
		  maskedHisto->setTitle("MASKED Pixel Index Hit Map;X Index [#];Y Index [#];Count [#]");
		}
        
		if(_doMaskedOutPlots){
		  streamlog_out ( MESSAGE4 ) <<  "unused channel on sensor "<<sensorID <<"... uMaxX: "<<uMaxX<<", uMaxY: "<<uMaxY<< std::endl;
		  tempHistoName = _maskedOutHistoName + "_d" + to_string( sensorID );
		  AIDA::IHistogram2D * maskedOutHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		  _maskedOutHistos.insert(std::make_pair(sensorID, maskedOutHisto));
		  maskedOutHisto->setTitle("MASKEDOUT Pixel Index Hit Map;X Index [#];Y Index [#];Count [#]");

		  if(_doLargeClusterMaps){
		    
		    tempHistoName = _maskedOutSelectedEventsHistoName + "_d" + to_string( sensorID );
		    //streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... xMax: "<<xMax<<", yMax: "<<yMax<< std::endl;
		    AIDA::IHistogram2D * maskedOutSelectedEventsHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		    _maskedOutSelectedEventsHistos.insert(std::make_pair(sensorID, maskedOutSelectedEventsHisto));
		    maskedOutSelectedEventsHisto->setTitle("MASKEDOUT Pixel Index Hit Map for Selected Events ;X Index [#];Y Index [#];Count [#]");
		  }
		}

		streamlog_out ( MESSAGE4 ) <<  "mapped sensor "<<sensorID <<"... maxX: "<<maxX<<", maxY: "<<maxY<< std::endl;
		tempHistoName = _mappedHistoName + "_d" + to_string( sensorID );
		xBin = maxX - minX + 1;
		xMin = static_cast<double >( minX ) - 0.5 ;
		xMax = static_cast<double >( maxX ) + 0.5;
		yBin = maxY - minY + 1;
		yMin = static_cast<double >( minY ) - 0.5;
		yMax = static_cast<double >( maxY ) + 0.5;
        //streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... xMax: "<<xMax<<", yMax: "<<yMax<< std::endl;
		AIDA::IHistogram2D * mappedHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		_mappedHistos.insert(std::make_pair(sensorID, mappedHisto));
		mappedHisto->setTitle("MAPPED Pixel Index Hit Map;X Index [#];Y Index [#];Count [#]");

		if(_doLargeClusterMaps){
		  tempHistoName = _mappedSelectedEventsHistoName + "_d" + to_string( sensorID );
		  //streamlog_out ( MESSAGE4 ) <<  "Sensor "<<sensorID <<"... xMax: "<<xMax<<", yMax: "<<yMax<< std::endl;
		  AIDA::IHistogram2D * mappedSelectedEventsHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), xBin, xMin, xMax,yBin, yMin, yMax);
		  _mappedSelectedEventsHistos.insert(std::make_pair(sensorID, mappedSelectedEventsHisto));
		  mappedSelectedEventsHisto->setTitle("MAPPED Pixel Index Hit Map for Selected Events;X Index [#];Y Index [#];Count [#]");
		}

		tempHistoName = _eventMultiplicityHistoName + "_d" + to_string( sensorID );
		int     eventMultiNBin  = 500;
		double  eventMultiMin   = 0;
		double  eventMultiMax   = 500;
		std::string  eventMultiTitle = "Event multiplicity";
		if ( isHistoManagerAvailable ) {
		  histoInfo = histoMgr->getHistogramInfo(  _eventMultiplicityHistoName );
		  if ( histoInfo ) {
		    streamlog_out ( DEBUG2 ) << (* histoInfo ) << std::endl;
		    eventMultiNBin  = histoInfo->_xBin;
		    eventMultiMin   = histoInfo->_xMin;
		    eventMultiMax   = histoInfo->_xMax;
		    if ( histoInfo->_title != "" ) eventMultiTitle = histoInfo->_title;
		  }
		}
		AIDA::IHistogram1D * eventMultiHisto =
		  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
		                                                            eventMultiNBin, eventMultiMin, eventMultiMax);
		_eventMultiplicityHistos.insert( std::make_pair(sensorID, eventMultiHisto) );
		eventMultiHisto->setTitle( eventMultiTitle.c_str() );

		if(_doMaskedOutPlots){// Spectrum of pixel charge of MASKEDOUT pixels
		  tempHistoName = _maskedOutSpectrumHistoName + "_d" + to_string( sensorID );
		  int    clusterNBin  = 61;
		  double clusterMin   = -0.5;
		  double clusterMax   = 60.5;
		  std::string clusterTitle = "Spectrum of pixel charge of MASKEDOUT pixels (in Detector specific Charge Unit);Charge;Count [#]";
		AIDA::IHistogram1D * maskedOutSpectrumHisto =
		  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
		                                                            clusterNBin,clusterMin,clusterMax);
		_maskedOutSpectrumHistos.insert( std::make_pair(sensorID, maskedOutSpectrumHisto) );
		maskedOutSpectrumHisto->setTitle( eventMultiTitle.c_str() );
		}

	}

  streamlog_out ( DEBUG5 )  << "end of Booking histograms " << std::endl; 
}
#endif
