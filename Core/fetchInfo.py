import re;

# The distance limit for what is considered a contact.
# This is based on this paper, which indicates that
# 7 angstroms is the best contact distance for
# distinguishing between folds:
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-292
DISTANCE_THRESHOLD = 7;
# Defining constants for fetchInfo returns.
i_PAR_NAME = 0;
i_PAR_AA = 1;
i_PAR_TYPE = 2;
i_PAR_RES = 3;
i_PAR_SS = 4;
i_CONTACT_NAME = 5;
i_CONTACT_AA = 6;
i_CONTACT_TYPE = 7;
i_CONTACT_RES = 8;
i_CONTACT_SS = 9;
i_CONTACT_AREA = 10;
i_SECTOR_NUM = 11;
i_RES_DIFF = 12;

# Tracking residue adjustments to enforce residues starting at 1.
resAdj = None;
adjusted = False;

# Purpose: To fetch info from a given line.
# Parameters:
#   line = The line containing the info.
#   ssList = A list of secondary structure values for the current residue.
#   i_ss = The anticipated index of the secondary structure element for this particle.
# Return:
#   A list of parameters fetched from the line.
def fetchLineInfo(line, ssList):
    global resAdj, adjusted;
    lineList = line.rstrip().lstrip().split(); # Splitting by any empty space.
    # Skipping non-informative lines. Also skipping solvent for now.
    if "#" in line or line == "\n" or line == "" or "SAS" in line or lineList[0] == "expected":
        return False;
    # Information about the source particle.
    parName = lineList[1];  # Name of the particle (V0, for example)
    parAA = parName[0];  # One-letter AA code for the particle
    parType = "BB" if parName[-1] == "0" else "SC";  # 0 indicates backbone, 1 indicates sidechain
    parRes = int(lineList[2]);  # Residue number of the particle (Like 1)
    if not adjusted:
        resAdj = parRes - 1;
        adjusted = True;
    parRes -= resAdj;
    parSS = ssList[parRes - 1];
    # Information about the contact particle.
    contactName = lineList[5];  # Name of the contact particle (L0, for example)
    contactAA = contactName[0];  # One-letter AA code for the contact particle
    contactType = "SV" if contactName == "O0" else "BB" if contactName[-1] == "0" else "SC"; # BB, SC, or SV.
    contactRes = int(lineList[6]) - resAdj if contactType != "SV" else "-"; # Residue of the contact particle (Like 2). 0 for SV.
    contactSS = ssList[contactRes - 1] if contactRes != "-" else "S";
    # General contact information.
    contactArea = float(lineList[8]);  # Area of the contact (>25 indicates a contact)
    sectorInfo = lineList[-1];  # Sector of the contact particle relative to the current particle
    # Calculating the difference between the residue positions.
    resDiff = abs(contactRes - parRes) if contactType != "SV" else "-";
    # Getting sector information.
    if "=" in sectorInfo:
        sectorNum = int(re.search('(?<=\=).+', sectorInfo).group(0).rstrip().lstrip());
        distanceInfo = lineList[-2];
    else:
        sectorNum = int(sectorInfo.rstrip().lstrip());
        distanceInfo = lineList[-3];
    if "=" in distanceInfo:
        contactDistance = float(re.search('(?<=\=).+', distanceInfo).group(0).rstrip().lstrip());
    else:
        contactDistance = float(distanceInfo);
    # Skipping lines describing sector 0 contacts or those with distances
    # below the threshold. Also skipping abnormal amino-acids.
    if sectorNum == 0 or contactDistance > DISTANCE_THRESHOLD or parAA == "X" or contactAA == "X":
        return False
    return (parName, parAA, parType, parRes, parSS,
            contactName, contactAA, contactType, contactRes, contactSS,
            contactArea, sectorNum, resDiff);

# Purpose: To filter a contact set, removing any contacts which
# share a sector with another, keeping the contact with the higher
# area.
# Parameters:
#   contactSet = The set of contacts for a particle.
# Return:
#   filterSet = The filtered set of contacts.
def filterSet(contactSet):
    toRemove = set();
    i = 0;
    while i < len(contactSet):
        contact = contactSet[i];
        otherContacts = [x for x in contactSet if x != contact];
        for otherContact in otherContacts:
            if contact[i_SECTOR_NUM] == otherContact[i_SECTOR_NUM]:
                toPop = contact if contact[i_CONTACT_AREA] < otherContact[i_CONTACT_AREA] else otherContact;
                toRemove.add(toPop);
        i += 1;
    filterSet = [x for x in contactSet if x not in toRemove];
    return filterSet;

# Purpose: To reset adjustment value and status.
def resetAdj():
    global resAdj, adjusted;
    resAdj = None;
    adjusted = False;