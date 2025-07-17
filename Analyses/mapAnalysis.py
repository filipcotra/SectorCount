import os;

# edge = (parName, parRes, parSector, parSS, contactName, contactRes, contactSector, contactSS);
i_PAR_NAME = 0;
i_PAR_RES = 1;
i_PAR_SECTOR = 2;
i_PAR_SS = 3;
i_CONTACT_NAME = 4;
i_CONTACT_RES = 5;
i_CONTACT_SECTOR = 6;
i_CONTACT_SS = 7;
# Defining constants for CASP rr maps.
LOWER_LIM = 0;
UPPER_LIM = 8;
PROB = 0.00;

# Purpose: To create contact maps from a given set of edges.
# Parameters:
#   edgeSet = The set of edges describing the structure.
#   aaSeq = The amino-acid sequence of the structure.
#   pdb = The pdb code for the structure.
#   outputDir = The directory where output will be stored.
def createMaps(edgeSet, aaSeq, pdb, outputDir):
    seqList = [];
    for i in aaSeq:
        seqList.append(aaSeq[i]);
    seqStr = "".join(seqList);
    outputName = f"{outputDir}/{pdb}_CM.txt";
    # Making headers.
    outputLine = f"##\t{seqStr}\n#\tExperimental Map - {pdb}\n";
    # Storing maps in dictionaries.
    BB2BB_map = set();
    SC2SC_map = set();
    SC2BB_map = set();
    ALL_map = set();
    # Iterating through the edges of the edge set, populating the dictionaries.
    for edge in edgeSet:
        # Extracting info from the sets.
        srcPar = edge[i_PAR_NAME];
        srcRes = int(edge[i_PAR_RES]) if edge[i_PAR_RES] != "-" else 0;
        dstPar = edge[i_CONTACT_NAME];
        dstRes = int(edge[i_CONTACT_RES]) if edge[i_CONTACT_RES] != "-" else 0;
        # Finding AA types.
        srcType = srcPar[-1];
        dstType = dstPar[-1];
        # Skipping solvent contacts. No self contacts. Keeping trivial contacts,
        # where separation is one or two.
        if srcPar == "O0" or dstPar == "O0" or srcRes == dstRes:
            continue;
        # Populating contact map dictionaries.
        contactKey = (srcRes, dstRes);
        revKey = (srcRes, dstRes);
        if srcType == "0" and dstType == "0": # BB2BB
            if revKey not in BB2BB_map:
                BB2BB_map.add(contactKey);
        elif srcType == "1" and dstType == "1": # SC2SC
            if revKey not in SC2SC_map:
                SC2SC_map.add(contactKey);
        else: #SC2BB
            if revKey not in SC2BB_map:
                SC2BB_map.add(contactKey)
        if revKey not in ALL_map:
            ALL_map.add(contactKey);
    # Outputting in CASP rr format.
    os.makedirs(f"{outputDir}/{pdb}", exist_ok = True);
    with open(f"{outputDir}/{pdb}/{pdb}.BB2BB.rr", "w") as file:
        file.write(buildFile_CASP(seqStr, BB2BB_map));
    with open(f"{outputDir}/{pdb}/{pdb}.SC2SC.rr", "w") as file:
        file.write(buildFile_CASP(seqStr, SC2SC_map));
    with open(f"{outputDir}/{pdb}/{pdb}.SC2BB.rr", "w") as file:
        file.write(buildFile_CASP(seqStr, SC2BB_map));
    with open(f"{outputDir}/{pdb}/{pdb}.rr", "w") as file:
        file.write(buildFile_CASP(seqStr, ALL_map));

# Purpose: To build a file in CASP rr format based on the provided
# contact information.
# Parameters:
#   aaSeq = The amino-acid sequence of the protein.
#   pairwiseProbs = A dictionary mapping residue position pairs to
#                   probability values.
# Return:
#   fileString = A formatted string representing the contact map.
def buildFile_CASP(aaSeq, pairwiseProbs):
    fileString = f"{aaSeq}\n";
    for pair in pairwiseProbs:
        resA, resB = pair;
        fileString += f"{int(resA)}\t{int(resB)}\t{LOWER_LIM}\t{UPPER_LIM}\t{PROB}\n";
    return fileString;

