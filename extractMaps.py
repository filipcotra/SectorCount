import sys;
import re;
from Analyses.mapAnalysis import createMaps;
from Core.fetchInfo import fetchLineInfo, filterSet, resetAdj;
from Core.extractEdgeSet import getEdgeSet;

# Defining constants for fetchInfo returns.
i_PAR_NAME = 0;
i_PAR_RES = 3;

# Only input should be a file with a list of files indicating
# contact results to turn into contact maps.
files = open(sys.argv[sys.argv.index("-i") + 1], "r");
ssDir = sys.argv[sys.argv.index("-s") + 1].rstrip("\\");
outputDir = sys.argv[sys.argv.index("-o") + 1] if "-o" in sys.argv else "";

# Looping through the file names.
for fileName in files:
    pdb = re.search(string = fileName, pattern = "(?<=contactResults\/).*(?=\.pdb)").group(0);
    aaSeq = {};
    # Getting information from the file.
    print(f"Creating Contact Map: {fileName}");
    currFile = open(fileName.rstrip(), "r");
    # Getting the equivalent ss file.
    pdbCode = re.search(".{4}(?=\.pdb)", fileName).group(0);
    ssFile = open(f"{ssDir}/{pdbCode}.ss", "r");
    ssLines = ssFile.readlines();
    ssInf = ssLines[-1].rstrip();
    ssList = list(ssInf);
    # Tracking the contact set of each particle for the current file.
    particleContacts = {};
    # In this file, identifying the sectors for non-backbone contacts and
    # adding the respective amino acids to the sector values. Fetching info
    # by particle.
    skipFile = False;
    lastLine = False;
    currPar = None;
    lastPar = None;
    contactSet = []; # Want to track all the contact info for the current particle.
    while not lastLine:
        line = currFile.readline();
        if "Warning" in line:
            skipFile = True;
            break;
        # Checking if we are at the end of the line.
        currentPos = currFile.tell(); # Identifying the current position of the file.
        lastLine = True if currFile.readline() == "" else False;
        currFile.seek(currentPos); # Resetting the file position to the current line.
        # If we run into any exceptions, which should be rare, just
        # continue to the next.
        try:
            lineInfo = fetchLineInfo(line, ssList);
        except:
            lineInfo = False;
        # If lineInfo is False, continue to the next line.
        # However, if we have reached the last line, don't.
        if not lineInfo and not lastLine:
            continue;
        # As long as lineInfo isn't false, we should collect new information.
        if lineInfo != False:
            # Importantly, this key uses the int version of the
            # particle residue number.
            currPar = (lineInfo[i_PAR_NAME], lineInfo[i_PAR_RES]);
            # Tracking the amino-acid sequence using the sequence number
            # as a key value and the residue one-letter code as the value.
            if currPar[1] not in aaSeq:
                aaSeq[currPar[1]] = currPar[0][0];
            if lastPar is None:
                lastPar = currPar;
            # If the current particle is different from the last, push the
            # current contact set.
            if currPar != lastPar:
                filteredSet = filterSet(contactSet);
                # Noting the contact set and then resetting it to empty.
                particleContacts[lastPar] = filteredSet;
                contactSet = [];
                lastPar = currPar;
            # Always adding the current line to the new contact set.
            contactSet.append(lineInfo);
        # If we are on the last line, push the current contact set.
        if lastLine:
            filteredSet = filterSet(contactSet);
            particleContacts[currPar] = filteredSet;
    currFile.close();
    if not skipFile:
        # Now we have all the contact sets for each particle
        # in the file, which we should extract to a set of overall
        # edges for the dataset.
        edgeSet = getEdgeSet(particleContacts);
        # Contact map analysis will create a contact map based on
        # the edgeSet and the aaSeq.
        createMaps(edgeSet, aaSeq, pdb, outputDir);
        # Deleting some of the bigger objects to hopefully save memory.
        del contactSet, particleContacts, edgeSet;
        resetAdj();
    else:
        print(f"{fileName} Skipped - Bad Data")
files.close();