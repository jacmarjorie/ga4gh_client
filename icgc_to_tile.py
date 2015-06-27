import sys
from threading import Thread, active_count 

# Append path to import file2tile 
sys.path.append('../tiledb')
from csvline import CSVLine
from file2tile import File2Tile, IDX, File2TileGlobals 
from bookkeeping import BookKeeping 

def progressPrint(outString):
  sys.stdout.write("\r" + outString)
  sys.stdout.flush()

class ICGC(File2Tile):
  """
  Class implementing the conversion from ICGC to Tile DB CSV
  """
  # reference to super object
  m_super = None

  # placeholder for values from previous line
  prev_SampleId = None
  prev_ChromosomePosition = None
  prev_TileDBPosition = None
  prev_TileDBValues = None
  prev_CallSetName = None
  prev_VariantSetId = None
  prev_VariantSetName = None
  GT = None

  # instance of the CSVLine object
  m_csv_line = None

  # Instance of BookKeeping
  m_BookKeeping = None

  def __init__(self):
    """
    Constructor takes the config file for the given input data
    """
    self.m_super = super(ICGC, self)
    self.m_super.__init__()

    self.m_csv_line = CSVLine()

    self.m_BookKeeping = ICGCGlobals.m_BookKeeping # BookKeeping()

  def generateCSV(self, inputFile, outFile, bVerbose = False):
    """
    Implements the ICGC specifics to generate a CSV output
    """
    self.initFilePointers(inputFile, outFile)

    nLine = 0L
    try:
      while( self.parseNextLine() ):
        self.checkSample()
        nLine += 1
        if( bVerbose ):
          progressPrint(str(nLine))

      # Cleanup and write the last CSV Line
      self.writeCSVLine()
    except:
      print "\nFAILED execution @ line " + str(nLine)
      raise
    finally: 
      # Close File Pointers
      self.closeFilePointers()

    # Reset to None in case the object is called to Generate CSV again
    self.prev_SampleId = None

  def checkSample(self):
    """
    Compares the current sample and the stored previous sample to return a bool
    - If a new sample is identified, write the existing data, and copy current
    - If it is the same sample, then compare with previous and append if new
    """
    self.VariantSetName = ICGCGlobals.variantMap\
                 [self.getValue("project_code")][ICGCGlobals.variantSetNameIdx]
    self.GT = self.TileDBValues["GT"]
    if( self.prev_SampleId == None ):
      # first sample copy data
      self.saveCurrentData()
    elif( self.prev_SampleId != self.SampleId or
          self.prev_CallSetName != self.CallSetName or
          self.prev_TileDBPosition != self.TileDBPosition or
          self.prev_VariantSetName != self.VariantSetName or 
          self.prev_TileDBValues["GT"] != self.GT ):
      # We have to write the current line and save a new line since the current
      # data differs by the SampleId, CallSetName, VariantSetName or TileDBPosition
      self.writeCSVLine()
      self.saveCurrentData()
    else:
      # This case is reached if the SampleId, CallSetName, TileDBPosition are 
      # the same. Verify if there is a delta in the fields whose type is an 
      # array, and append new data if it is not already in the prev_TileDBValues
      for key in self.TileDBValues.keys():
        if( key in CSVLine.arrayFields ):
          for value in self.TileDBValues[key]:
            if( value not in self.prev_TileDBValues[key] ):
              self.prev_TileDBValues[key].append(value)

  def saveCurrentData(self):
    """
    Copies the current data into the prev_ variables
    """
    self.prev_SampleId = self.SampleId
    self.prev_CallSetName = self.CallSetName
    self.prev_ChromosomePosition = self.ChromosomePosition
    self.prev_TileDBPosition = self.TileDBPosition
    self.prev_TileDBValues = self.TileDBValues
    # If the Variant Set Name changes then update the Varient Set Id
    # VariantSetName changes very infrequently so this optimization will help 
    if( self.prev_VariantSetName != self.VariantSetName ):
      self.prev_VariantSetId = self.m_BookKeeping.getVariantSetId(self.VariantSetName)
      self.prev_VariantSetName = self.VariantSetName

  def writeCSVLine(self):
    """
    Creates the CSV object from the data set we have and writes it to disk
    """
    if( self.prev_SampleId == None ):
      return
    """
    Since ICGC represents insertions and deletions with a '-' it needs to 
    be processed before adding to the CSV Line
    1. Insertions
      Ref = get_ref_from_ensembl(ICGC_chr, ICGC_start-1, ICGC_start - 1)
      Alt = Ref + Alt # Concatenate Ref and Alt
      Chr = ICGC_Chr
      Start = ICGC_Start - 1
      End = Start 
    2. Deletions
      NewRef = get_ref_from_ensembl(ICGC_chr, ICGC_start-1, ICGC_start-1) 
      Ref = NewRef + ICGC_Ref
      Alt = NewRef
      Chr = ICGC_Chr
      Start = ICGC_Start - 1
      End = ICGC_End
    """
    # Insertion
    if( self.prev_TileDBValues["REF"] == "-" ):
      assembly = self.prev_ChromosomePosition[IDX.CHR_ASSEMBLY]
      chromosome = self.prev_ChromosomePosition[IDX.CHR_CHR]
      start = self.prev_ChromosomePosition[IDX.CHR_START] - 1
      end = start # self.prev_ChromosomePosition[IDX.CHR_END]

      ref = self.translator.getReference(assembly, chromosome, start, start)
      self.prev_TileDBValues["REF"] = ref 
      index = 0
      for value in self.prev_TileDBValues["ALT"]:
        self.prev_TileDBValues["ALT"][index] = ref + value 
        index += 1
      self.prev_TileDBPosition = self.translator.getTileDBPostition(assembly, 
                                                      chromosome, start, end)
    else:
      bFlag = False
      for value in self.prev_TileDBValues["ALT"]:
        if( value == "-" ):
          bFlag = True
          break
      if( bFlag ):
        assembly = self.prev_ChromosomePosition[IDX.CHR_ASSEMBLY]
        chromosome = self.prev_ChromosomePosition[IDX.CHR_CHR]
        start = self.prev_ChromosomePosition[IDX.CHR_START] - 1
        end = self.prev_ChromosomePosition[IDX.CHR_END]

        ref = self.translator.getReference(assembly, chromosome, start, start)
        self.prev_TileDBValues["REF"] = ref + self.prev_TileDBValues["REF"]
        index = 0
        for value in self.prev_TileDBValues["ALT"]:
          if( value == "-" ):
            self.prev_TileDBValues["ALT"][index] = ref
          else:
            self.prev_TileDBValues["ALT"][index] = ref + value 
          index += 1

        self.prev_TileDBPosition = self.translator.getTileDBPostition(assembly, 
                                                        chromosome, start, end)

    # get the call set id which is also the Sample Id or Row Id for Tile DB
    prev_CallSetId = self.m_BookKeeping.getCallSetId(self.prev_SampleId, 
                                                self.prev_CallSetName)
    # Update the book keeping structure with the VariantSetId to CallSetId map
    self.m_BookKeeping.updateVariantSetId(self.prev_VariantSetId, prev_CallSetId)

    self.m_csv_line.set("SampleId", prev_CallSetId) 
    self.m_csv_line.set("Location", self.prev_TileDBPosition[IDX.START])
    self.m_csv_line.set("End", self.prev_TileDBPosition[IDX.END])

    # set ALT value first to avoid failures from set function
    # Append "&" at the very end since that is what TileDB expects
    self.m_csv_line.set("ALT", self.prev_TileDBValues["ALT"])
    self.m_csv_line.set("PLOIDY", ICGCGlobals.cPLOIDY)

    for key in self.prev_TileDBValues.keys():
      if( key == "ALT" ):
        continue
      value = self.prev_TileDBValues[key]

      # For ICGC the GT value from the MAF file is the first entry and
      # the second is always a 0
      if( key == "GT" ):
        if( value[0] != '1' ):
          print value[0]
        value.append('0')
      # Cleanup Value
      if( value == "" ):
        value = "#"
      self.m_csv_line.set(key, value)
    # Write the line into file
    self.outFile.write(self.m_csv_line.getCSVLine() +  "\n")

class ICGCGlobals:
  m_BookKeeping = BookKeeping()

  variantMap = dict()
  variantSetNameIdx = 2

  # For ICGC the PLOIDY value is always 2
  cPLOIDY = 2

  def initICGCGlobals(self, variantConfig):
    """
    Takes an input config file that contains the ICGC variant mapping
    """
    import json
    fp = open(variantConfig, 'r')
    ICGCGlobals.variantMap = json.load(fp)
    fp.close()

def getFileName(inFile, splitStr = None):
  """
  Strips the /'s and gets the file name without extension
  """
  paths = inFile.strip().split('/')

  if( splitStr == None ):
    splitStr = "."
  fileNames = paths[len(paths) - 1].split(splitStr)
  return ".".join(fileNames[0:len(fileNames) - 1])

def poolGenerateCSV(config, inputFile, outFile):
  """
  This function is used by multiprocess.Pool to generate CSV for each input file
  """
  ic = ICGC()
  ic.generateCSV(inputFile, outFile)

def getFilePointer(fileName, type, mode):
  """
  returns a file pointer given the name, type and mode 
  """
  if( type == "tsv" ):
    return open(fileName, mode)
  elif( type == "tsv.gz" ):
    import gzip
    return gzip.open(fileName, mode)
  else:
    return open(fileName, mode)

def parallelGen(config, inputFileList, outputDir, combinedOutputFile, fileType):
  """
  Function that spawns the Pool of ICGC objects to work on each of the input files
  Once the BookKeeping support moves to a real DB, move from threads to multiprocessing.Pool
  """
  m_threads = list()

  count = 0

  jsonMap = File2TileGlobals.initFile2TileGlobals(File2TileGlobals(), config)
  ICGCGlobals.initICGCGlobals(ICGCGlobals(), jsonMap["VariantConfig"])

  for inFile in inputFileList:
    count += 1
    outFile = outputDir + "/" + getFileName(inFile, "." + fileType) + ".csv"
    m_thread = Thread(target=poolGenerateCSV, 
                      name="Thread " + str(count) + " :" + inFile, 
                      args=((config, getFilePointer(inFile, fileType, "r"), 
                            getFilePointer(outFile, "csv", "w"))))
    m_thread.start()
    m_threads.append((m_thread, outFile))

  import os
  combinedOutput = outputDir + "/" + combinedOutputFile
  countCombined = 0
  isCombined = [False] * len(m_threads)
  progressPrint("%d of %d Completed"%(countCombined, len(m_threads))) 
  while( active_count() > 0 and countCombined != len(m_threads) ) :
    index = -1
    for (m_thread, outFile) in m_threads:
      index += 1
      if( isCombined[index] ):
        continue
      if( m_thread.is_alive() ):
        m_thread.join(120)
      else:
        os.system('cat %s >> %s'%(outFile, combinedOutput))
        countCombined += 1
        isCombined[index] = True
        progressPrint("%d of %d Completed"%(countCombined, len(m_threads))) 
  print ""

  jsonOut = outputDir + "/" + getFileName(combinedOutputFile) + ".json"
  fp = open(jsonOut, "w")
  fp.write(ICGCGlobals.m_BookKeeping.getJSON())
  fp.close()

if __name__ == "__main__":
  import argparse 
  parser = argparse.ArgumentParser(description = "Convert ICGC MAF format to Tile DB CSV")

  parser.add_argument("-c", "--config", required=True, type=str, 
                      help="input configuration file for ICGC conversion")
  parser.add_argument("-o", "--output", required=True, type=str, 
                      help="output Tile DB CSV file (without the path) which will be stored in the output directory")
  parser.add_argument("-d", "--outputdir", required=True, type=str, 
                      help="Output directory where the outputs need to be stored")
  parser.add_argument("-i", "--inputs", nargs='+', type=str, required=True, 
                      help="List of input ICGC files to convert")
  parser.add_argument("-t", "--inputtype", required=True, choices=['tsv', 'tsv.gz'], 
                      help="Type of input file")
  parser.add_argument("-a", "--appendJSON", required=False, default = "None", 
                      help="JSON file to which the new data should be append. The output JSON will be the new appended JSON, and the original JSON will be left intact.")

  args = parser.parse_args()

  if( args.appendJSON != "None" ):
    ICGCGlobals.m_BookKeeping.loadJSON(args.appendJSON)

  parallelGen(args.config, args.inputs, args.outputdir, args.output, args.inputtype)
