from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport numpy as np

cdef extern from "VBF/Words.h":
  pass
# From Words.h
ctypedef np.int8_t byte
ctypedef np.int16_t word16
ctypedef np.int32_t word32
ctypedef np.int64_t word64
ctypedef np.uint8_t ubyte
ctypedef np.uint16_t uword16
ctypedef np.uint32_t uword32
ctypedef np.uint64_t uword64      

# Cython doesn't support embedded enums in c++ classes yet,
# so we'll not use the cppVEventType, and use a pure python
# extension class instead.  This is fine for read only.


cdef extern from "VBF/VDatum.h":
  cdef cppclass cppVEvent "VEvent":
    cppVEvent()
    cppVEvent* copyEvent()
    uword32 getSize()
    uword16 getNumSamples()
    uword16 getNumChannels()
    uword16 getMaxNumChannels()
    uword16 getNumClockTrigBoards()
    bool getCompressedBit()
    void resizeChannelData(uword16, uword16)
    uword32* getHitPattern()
    bool getHitBit(uword32)
    uword32* getTriggerPattern()
    bool getTriggerBit(unsigned)
    uword16 getPedestalAndHiLo(unsigned)
    uword16 getPedestalAndHiLoWithoutVerify(unsigned)
    uword16 getPedestal(unsigned)
    bool getHiLo(unsigned)
    uword16 getCharge(unsigned)
    uword16 getChargeWithoutVerify(unsigned)
    ubyte* getSamplePtr(unsigned, unsigned)
    ubyte getSample(unsigned, unsigned)
    uword32* getClockTrigData(unsigned)
    vector[bool] getFullHitVec()
    vector[bool] getFullTrigVec()
    # from VDatum
    unsigned getGPSTimeNumElements()
    ubyte getNodeNumber()
    ubyte getTriggerMask()
    uword32 getEventNumber()
    uword16* getGPSTime()
    ubyte getGPSYear()
    ubyte getRawEventTypeCode()
    ubyte getEventTypeCode()
  cdef cppclass cppVArrayTrigger "VArrayTrigger":
    cppVArrayTrigger()
    uword32 getSize()
    ubyte getNumSubarrayTelescopes()
    ubyte getNumTriggerTelescopes()
    uword16 getATFlags()
    ubyte getConfigMask()
    uword32 getRunNumber()
    bool hasTenMHzClockArray()
    bool hasOptCalCountArray() 
    bool hasPedCountArray()
    uword32* getTenMHzClockArray()
    uword32* getOptCalCountArray()
    uword32* getPedCountArray()
    uword32 getSubarrayTelescopeId(unsigned)
    uword32 getTDCTime(unsigned)
    uword32 getSpecificRawEventTypeCode(unsigned)
    uword32 getSpecificEventTypeCode(unsigned)
    float getAzimuth(unsigned)
    float getAltitude(unsigned)
    uword32 getShowerDelay(unsigned)
    bool hasShowerDelay()
    uword32 getCompDelay(unsigned)
    bool hasCompDelay()
    bool hasL2CountsArray()
    bool hasCalCountsArray()
    uword32* getL2CountsArray(unsigned)
    uword32* getCalCountsArray(unsigned)
    uword32 getDelay(unsigned)
    bool hasDelay()
    uword32 getL2ScalarRate(unsigned)
    uword32* getL2Pattern(unsigned)
    bool hasL2Pattern()
    uword32 getTenMhzClock(unsigned)
    bool hasTenMhzClock()
    uword32 getVetoedClock(unsigned)
    uword32 getTriggerTelescopeId(unsigned)
    # from VDatum
    unsigned getGPSTimeNumElements()
    ubyte getNodeNumber()
    ubyte getTriggerMask()
    uword32 getEventNumber()
    uword16* getGPSTime()
    ubyte getGPSYear()
    ubyte getRawEventTypeCode()
    ubyte getEventTypeCode()

cdef extern from "VBF/VArrayEvent.h":
  cdef cppclass cppVArrayEvent "VArrayEvent":
    cppVArrayEvent() except +
    unsigned getNumEvents()
    unsigned getMaxNumEvents()
    unsigned getNumDatums()
    long getRun()
    bool hasEventNumber()
    unsigned long getEventNumber() except +
    bool hasTrigger()
    cppVArrayTrigger* getTrigger() except +
    cppVEvent* getEventAt(unsigned) except +
    cppVEvent* getEvent(unsigned) except +
    cppVEvent* getEventByNodeNumber(unsigned) except +
    bool isEmpty()
    vector[bool] getExpectedTelescopes()
    vector[bool] getPresentTelescopes()
    string getSummary() except +

cdef extern from "VBF/VPacket.h":
  cdef cppclass cppVPacket "VPacket":
    cppVPacket()
    bool empty()
    unsigned size()
    bool hasArrayEvent()
    bool hasSimulationHeader()
    bool hasSimulationData()
    bool hasEventOverflow()
    cppVArrayEvent* getArrayEvent() except +

cdef extern from "VBF/VBankFileReader.h":
  cdef cppclass cppVBankFileReader "VBankFileReader":
    cppVBankFileReader(const string &filename, bool map_index, bool read_only) except +
    long getRunNumber()
    bool hasIndex()
    uword32 numPackets()
    bool hasPacket(uword32)
    cppVPacket* readPacket(uword32) except +
    uword32 generateIndexAndChecksum()
    uword32 getChecksum()
    uword32 calculateChecksum()
    uword64 getFileSize()


