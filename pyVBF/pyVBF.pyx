from cython.operator cimport dereference as deref
from cython cimport view
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport numpy as np
import numpy as np

# Cython doesn't support embedded enums in c++ classes yet,
# so we'll not use the cppVEventType, and use a pure python
# extension class instead.  This is fine for read only.

cdef class VEventType:
  """Python workaline of VBF VEventType

  Like the C++ version, this class is intialized to default values
  and the properties can be set individually or by setNewStyleCode()
  """
  cdef readonly object TriggerType
  cdef readonly object CalibrationType      
  cdef int trigger, calibration
  cdef bool xtra_samples, force_full_mode, has_error
  def __cinit__(self):
    self.TriggerType = ('L2_TRIGGER', 'HIGH_MULT_TRIGGER', 'NEW_PHYS_TRIGGER', 'CAL_TRIGGER', 'PED_TRIGGER')
    self.CalibrationType = ('NOT_CALIBRATION', 'OPTICAL_CALIBRATION', 'CHARGE_CALIBRATION')      
  def __init__(self):
    self.trigger = 0
    self.xtra_samples = False
    self.force_full_mode = False
    self.calibration = 0
    self.has_error = False
  property trigger:  
    def __get__(self):
      return self.TriggerType[self.trigger]
    def __set__(self, int value):
      try:
        self.trigger = self.TriggerType.index(value)
      except:
        self.trigger = value
  property calibration:
    def __get__(self):
      return self.CalibrationType[self.calibration]
    def __set__(self, int value):
      try:
        self.calibration = self.CalibrationType.index(value)
      except:
        self.calibration = value
  property xtra_samples:
    def __get__(self):
      return self.xtra_samples
    def __set__(self, bool value):
      self.xtra_samples = value
  property force_full_mode:
    def __get__(self):
      return self.force_full_mode
    def __set__(self, bool value):
      self.force_full_mode = value
  property has_error:
    def __get__(self):
      return self.has_error
  def setNewStyleCode(self, int code):
    """Set the event type according to a numerical code"""
    self.trigger = 0
    self.xtra_samples = False
    self.force_full_mode = False
    self.calibration = 0
    self.has_error = False
    if code<1 or code>14:
      self.has_error = True
      return
    if code in (4,5,6):
      self.trigger = 1
    elif code in (7,8,9):
      self.trigger = 2
    elif code == 10:
      self.trigger = 3
    elif code in (11,13):
      self.trigger = 4      
    else:
      pass
    if code in (3,6,9):
      self.force_full_mode = True
    elif code in (2,5,8):
      self.xtra_samples = True
    if code in (11,12):
      self.calibration = 1
    elif code in (13,14):
      self.calibration = 2      
    else:
      pass
  def getBestNewStyleCode(self):
    """Return a neumerical code corresponding to the set event type"""
    if self.trigger == 0:
      if self.calibration == 1:
        return 12
      elif self.calibration == 2:
        return 14
      else:
        if self.xtra_samples:
          if self.force_full_mode:
            return 3
          else:
            return 2
        else:
          return 1
    elif self.trigger == 1:
      if self.xtra_samples:
        if self.force_full_mode:
          return 6
        else:
          return 5
      else:
        return 4
    elif self.trigger == 2:
      if self.xtra_samples:
        if self.force_full_mode:
          return 9
        else:
          return 8
      else:
        return 7
    elif self.trigger == 3:
      if self.calibration == 2:
        return 13
      else:
        return 11
    elif self.trigger == 4:
      return 10
    else:
      return -1      

cdef class VEvent:
  """Trace information for an individual telescope"""
  cdef cppVEvent* __thisptr
  cdef int __delete
  cdef int __i # iterator index
  def __cinit__(self, *args):
    if len(args) == 2:
      self.__thisptr = new cppVEvent()
      self.__delete = 1
      self.__i = 0
  cdef __initFromRawPointer(self, cppVEvent* cptr):
    print(<long>cptr)
    self.__thisptr = cptr
    self.__delete = 0
  def __dealloc__(self):
    if self.__delete:
      del self.__thisptr
  def getNumSamples(self):
    """Number of samples per channel per event"""
    return self.__thisptr.getNumSamples()
  def getNumChannels(self):
    """Number of channels in this telescope"""
    return self.__thisptr.getNumChannels()
  def getMaxNumChannels(self):
    return self.__thisptr.getMaxNumChannels()
  def getNumClockTrigBoards(self):
    return self.__thisptr.getNumClockTrigBoards()
  def getCompressedBit(self):
    return self.__thisptr.getCompressedBit()
  def getHitPattern(self):
    return deref(self.__thisptr.getHitPattern())
  def getHitBit(self, int i):
    return self.__thisptr.getHitBit(i)
  def getTriggerPattern(self):
    return deref(self.__thisptr.getTriggerPattern())
  def getTriggerBit(self, int i):
    return self.__thisptr.getTriggerBit(i)
  def getPedestalAndHiLo(self, int i):
    return self.__thisptr.getPedestalAndHiLo(i)
  def getPedestalAndHiLoWithoutVerify(self, int i):
    return self.__thisptr.getPedestalAndHiLoWithoutVerify(i)
  def getPedestal(self, int i):
    return self.__thisptr.getPedestal(i)
  def getHiLo(self, int i):
    return self.__thisptr.getHiLo(i)
  def getCharge(self, int i):
    return self.__thisptr.getCharge(i)
  def getSample(self, int i, int j):
    return self.__thisptr.getSample(i,j)
  def getClockTrigData(self,i):
    return deref(self.__thisptr.getClockTrigData(i))
  # From VDatum
  def getNodeNumber(self):
    """Return the telescope number idexed from zero (2->T3)"""
    return self.__thisptr.getNodeNumber()
  def getTriggerMask(self):
    return self.__thisptr.getTriggerMask()
  def getEventNumber(self):
    return self.__thisptr.getEventNumber()
  def getGPSTime(self):
    cdef uword16* start
    cdef int i
    start = self.__thisptr.getGPSTime()
    return [deref(start+i) for i in range(self.__thisptr.getGPSTimeNumElements())]
    #return deref(self.__thisptr.getGPSTime())
  def getGPSYear(self):
    return self.__thisptr.getGPSYear()
  def getEventType(self):
    et = VEventType()
    et.setNewStyleCode(self.__thisptr.getEventTypeCode())
    return et
  def getRawEventTypeCode(self):
    return self.__thisptr.getRawEventTypeCode()
  def getEventTypeCode(self):
    return self.__thisptr.getEventTypeCode()
  def slowArray(self):
    """Return a [c,s] dimension arrary

    c is the number of channels,
    s is the number of samples

    NOTE: this function uses the public getSample() method of VEvent, and
    should always be expected to work.
    """
    cdef int c,s, nc, ns
    nc = self.getNumChannels()
    ns = self.getNumSamples()
    cdef np.ndarray[np.uint8_t, ndim=2] a = np.zeros([nc, ns], dtype=np.uint8)
    for c in range(nc):
      for s in range(ns):
        a[c,s] = self.getSample(c,s)
    return a
  def array(self,channel=None):
    """Return an arrary of samples.

    Keyword arguments:
    channel -- the requested channel (starts at 0) (default None)

    If the requested channel is None, all channels will be read out in 
    a [c,s] two dimensional array.

    c is the number of channels,
    s is the number of samples

    NOTE: this function uses pointers returned by public methods of VEvent
    but assumes the memory for all channels is continious.
    """
    cdef int nc, ns
    nc = self.getNumChannels()
    ns = self.getNumSamples()
    if channel is not None:
      return np.asarray(<np.uint8_t[:ns]> self.__thisptr.getSamplePtr(channel,0)).copy()
    else:
      return np.asarray(<np.uint8_t[:nc, :ns]> self.__thisptr.getSamplePtr(0,0)).copy()      
    # cdef ubyte[:,::1] view = <ubyte[:nc, :ns]> self.__thisptr.getSamplePtr(0,0)
    # return np.array(view)
  def HiLoArray(self):
    """Return a boolean array for the high-low bit for each channel

    False -> High gain
    True -> Low gain
    """
    cdef int c, nc
    nc = self.getNumChannels()
    # cython seems to not have a np.bool_t, so fake it.
    cdef np.ndarray[np.uint8_t, ndim=1] a = np.zeros([nc], dtype=np.uint8)
    for c in range(nc):
      a[c] = self.getHiLo(c)
    return a.astype(np.bool)
  # Special methods for python
  def __len__(self):
    return self.getNumChannels()   
  def __iter__(self):
    return self.fastArray()
  def __getitem__(self, int i):
    if i < self.getNumChannels():
      return self.array(i)
    else:
      raise IndexError   

cdef class VArrayTrigger:
  cdef cppVArrayTrigger* __thisptr
  cdef int __delete
  cdef int __i  # iterator index
  def __cinit__(self, *args):
    if len(args) == 2:
      self.__thisptr = new cppVArrayTrigger()
      self.__delete = 1
  cdef __initFromRawPointer(self, cppVArrayTrigger* cptr):
    self.__thisptr = cptr
    self.__delete = 0
  def __dealloc__(self):
    if self.__delete:
      del self.__thisptr
  def getNumSubarrayTelescopes(self):
    return self.__thisptr.getNumSubarrayTelescopes()
  def getNumTriggerTelescopes(self):
    return self.__thisptr.getNumTriggerTelescopes()
  def getATFlags(self):
    return self.__thisptr.getATFlags()
  def getConfigMask(self):
    return self.__thisptr.getConfigMask()
  def getRunNumber(self):
    return self.__thisptr.getRunNumber()
  def hasTenMHzClockArray(self):
    return self.__thisptr.hasTenMHzClockArray()
  def hasOptCalCountArray(self):
    return self.__thisptr.hasOptCalCountArray()
  def hasPedCountArray(self):
    return self.__thisptr.hasPedCountArray()
  def getTenMHzClockArray(self):
    return deref(self.__thisptr.getTenMHzClockArray())
  def getOptCalCountArray(self):
    return deref(self.__thisptr.getOptCalCountArray())
  def getPedCountArray(self):
    return deref(self.__thisptr.getPedCountArray())
  def getSubarrayTelescopeId(self, int i):
    return self.__thisptr.getSubarrayTelescopeId(i)
  def getTDCTime(self, int i):
    return self.__thisptr.getTDCTime(i)
  def getAzimuth(self, int i):
    return self.__thisptr.getAzimuth(i)
  def getAltitude(self, int i):
    return self.__thisptr.getAltitude(i)
  def getShowerDelay(self, int i):
    return self.__thisptr.getShowerDelay(i)
  def hasShowerDelay(self):
    return self.__thisptr.hasShowerDelay()
  def getCompDelay(self, int i):
    return self.__thisptr.getCompDelay(i)
  def hasCompDelay(self):
    return self.__thisptr.hasCompDelay()
  def hasL2CountsArray(self):
    return self.__thisptr.hasL2CountsArray()
  def hasCalCountsArray(self):
    return self.__thisptr.hasCalCountsArray()
  def getL2CountsArray(self, int i):
    return deref(self.__thisptr.getL2CountsArray(i))
  def getCalCountsArray(self, int i):
    return deref(self.__thisptr.getCalCountsArray(i))
  def getDelay(self, int i):
    return self.__thisptr.getDelay(i)
  def hasDelay(self):
    return self.__thisptr.hasDelay()
  def getL2ScalarRate(self, int i):
    return self.__thisptr.getL2ScalarRate(i)
  def getL2Pattern(self, int i):
    return deref(self.__thisptr.getL2Pattern(i))
  def hasL2Pattern(self):
    return self.__thisptr.hasL2Pattern()
  def getTenMhzClock(self, int i):
    return self.__thisptr.getTenMhzClock(i)
  def hasTenMhzClock(self):
    return self.__thisptr.hasTenMhzClock()
  def getVetoedClock(self, int i):
    return self.__thisptr.getVetoedClock(i)
  def getTriggerTelescopeId(self, int i):
    return self.__thisptr.getTriggerTelescopeId(i)
  # From VDatum
  def getNodeNumber(self):
    return self.__thisptr.getNodeNumber()
  def getTriggerMask(self):
    return self.__thisptr.getTriggerMask()
  def getEventNumber(self):
    return self.__thisptr.getEventNumber()
  def getGPSTime(self):
    return deref(self.__thisptr.getGPSTime())
  def getGPSYear(self):
    return self.__thisptr.getGPSYear()
  def getEventType(self):
    et = VEventType()
    et.setNewStyleCode(self.__thisptr.getEventTypeCode())
    return et
  def getRawEventTypeCode(self):
    return self.__thisptr.getRawEventTypeCode()
  def getEventTypeCode(self):
    return self.__thisptr.getEventTypeCode()      

cdef class VArrayEvent:      
  cdef cppVArrayEvent* __thisptr
  cdef int __delete
  cdef int __event # Index for iterator
  def __cinit__(self, *args):
    if len(args) == 2:
      self.__thisptr = new cppVArrayEvent()
      self.__delete = 1
      self.__event = 0
  cdef __initFromRawPointer(self, cppVArrayEvent* cptr):
    self.__thisptr = cptr
    self.__delete = 0
    self.__event = 0
  def __dealloc__(self):
    if self.__delete:
      del self.__thisptr
  def getNumEvents(self):
    return self.__thisptr.getNumEvents()
  def getMaxNumEvents(self):
    return self.__thisptr.getMaxNumEvents()
  def getRun(self):
    return self.__thisptr.getRun()
  def hasEventNumber(self):
    return self.__thisptr.hasEventNumber()
  def getEventNumber(self):
    return self.__thisptr.getEventNumber()
  def hasTrigger(self):
    return self.__thisptr.hasTrigger()
  def getTrigger(self):
    if not self.hasTrigger():
      raise Exception
    t = VArrayTrigger()
    cdef cppVArrayTrigger* ptr = self.__thisptr.getTrigger()
    t.__initFromRawPointer(ptr)
    return t
  def getEventAt(self, int i):
    if i > self.getNumEvents():
      raise IndexError
    e = VEvent()
    cdef cppVEvent* ptr = self.__thisptr.getEvent(i)
    e.__initFromRawPointer(ptr)
    return e
  def getEvent(self, int i):
    return self.getEventAt(i)
  def getEventByNodeNumber(self, int i):
    """Get event by telescope number indexed from zero (2->T3)"""
    e = VEvent()
    cdef cppVEvent* ptr = self.__thisptr.getEventByNodeNumber(i)
    if ptr == NULL:
      raise IndexError
    e.__initFromRawPointer(ptr)
    return e
  def isEmpty(self):
    return self.__thisptr.isEmpty()
  def getExpectedTelescopes(self):
    return self.__thisptr.getExpectedTelescopes()
  def getPresentTelescopes(self):
    return self.__thisptr.getPresentTelescopes()
  def getSummary(self):
    return self.__thisptr.getSummary()
  # Special Methods for python
  def __len__(self):
    return self.getNumEvents()
  def __iter__(self):
    return self
  def __next__(self):
    if self.__event < self.getNumEvents():
      event = self.getEventAt(self.__event)
      self.__event += 1
      return event
    else:
      raise StopIteration
  def __getitem__(self, int i):
    if i < self.getNumEvents():
      return self.getEventAt(i)
    else:
      raise IndexError      
  def __repr__(self):
    return "<"+str(self)+">"
  def __str__(self):
    return "VArrayEvent containting {0} telescope events".format(self.getNumEvents())

cdef class VPacket:      
  cdef cppVPacket* __thisptr
  cdef int __delete
  def __cinit__(self, *args):
    if len(args) == 2:
      self.__thisptr = new cppVPacket()
      self.__delete = 1
  cdef __initFromRawPointer(self, cppVPacket* cptr):
    self.__thisptr = cptr
    self.__delete = 0
  def __dealloc__(self):
    if self.__delete:
      del self.__thisptr
  property size:
    def __get__(self):
      return self.__thisptr.size()
  property empty:
    def __get__(self):
      return self.__thisptr.empty()
  def hasArrayEvent(self):
    return self.__thisptr.hasArrayEvent()
  def hasSimulationHeader(self):
    return self.__thisptr.hasSimulationHeader() 
  def hasSimulationData(self):
    return self.__thisptr.hasSimulationData()
  def hasEventOverflow(self):
    return self.__thisptr.hasEventOverflow()
  def getArrayEvent(self):
    vae = VArrayEvent(self)
    cdef cppVArrayEvent* ptr = self.__thisptr.getArrayEvent()
    vae.__initFromRawPointer(ptr)
    return vae
  def keys(self):
    keys = []
    if self.__thisptr.hasArrayEvent():
      keys.append("ArrayEvent")
    if self.__thisptr.hasSimulationHeader():
      keys.append("SimulationHeader")
    if self.__thisptr.hasSimulationData():
      keys.append("SimulationData")
    if self.__thisptr.hasEventOverflow():
      keys.append("EventOverflow")
    return keys
  def items(self):
    items = []
    for key in self.keys():
      items.append((key,self[key]))
    return items
  def values(self):
    values = []
    for key in self.keys():
      values.append(self[key])
    return values
  def __len__(self):
    return self.__thisptr.size()
  def __getitem__(self, key):
    if key == "ArrayEvent":
      return self.getArrayEvent()
    elif key == "SimulationHeader":
      print("not supported")
    elif key == "SimulationData":
      print("not supported")
    elif key == "EventOverflow":
      print("not supported")
    else:
      raise KeyError
  def __contains__(self,key):
    return key in self.keys()
  def __repr__(self):
    return "<"+str(self)+">"
  def __str__(self):
    if self.empty:
      return "Empty VPacket"
    return "VPacket with keys " + str(self.keys())

cdef class VBankFileReader:
  """Opens a .vbf (or .cvbf) file for reading."""
  cdef cppVBankFileReader* __thisptr  # hold a C++ instance which we're wrapping
  cdef int __packet
  def __cinit__(self, string filename, bool map_index=True, bool read_only=True):
    self.__thisptr = new cppVBankFileReader(filename,map_index,read_only)
    self.__packet = 0
  def __dealloc__(self):
    del self.__thisptr
  def getRunNumber(self):
    """Return the run number for this file.  

    A VBF file may contain data for only one run."""
    return self.__thisptr.getRunNumber()
  def hasIndex(self):
    return self.__thisptr.hasIndex()
  def numPackets(self):
    return self.__thisptr.numPackets()
  def readPacket(self, i=None):
    """Return a VPacket

    There is a VPacket for each event number from 0 to numPackets()-1
    
    If no packet number is specified, the next available packet will be
    returned.
    """
    if i is not None and self.hasIndex():
      self.__packet = i
    vpkt = VPacket()
    cdef cppVPacket* ptr = self.__thisptr.readPacket(self.__packet)
    self.__packet += 1
    vpkt.__initFromRawPointer(ptr)
    return vpkt
  # Special Methods for python
  def __len__(self):
    if self.hasIndex():
      return self.__thisptr.numPackets()
  def __iter__(self):
    return self
  def __next__(self):
    if self.__packet < self.numPackets():
      return self.readPacket()
    else:
      raise StopIteration
  def __getitem__(self, int i):
    if i < self.numPackets:
      return self.readPacket(i)
    else:
      raise IndexError
  def __repr__(self):
    return "<"+str(self)+">"
  def __str__(self):
    return "VBF File.  Run number: {0:d}, packets: {1:d}".format(self.getRunNumber(),self.numPackets()) 


