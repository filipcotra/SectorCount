import pandas as pd;
import random as rand;
from Core.fetchInfo import i_PAR_NAME, i_PAR_AA, i_PAR_TYPE, i_PAR_RES, i_PAR_SS,\
        i_CONTACT_NAME, i_CONTACT_AA, i_CONTACT_TYPE, i_CONTACT_RES, i_CONTACT_SS,\
        i_CONTACT_AREA, i_SECTOR_NUM, i_RES_DIFF;

# Defining constants for navigating edges.
i_EDGE_NAME_A = 0;
i_EDGE_RES_A = 1;
i_EDGE_SECTOR_A = 2;
i_EDGE_SS_A = 3;
i_EDGE_NAME_B = 4;
i_EDGE_RES_B = 5;
i_EDGE_SECTOR_B = 6;
i_EDGE_SS_B = 7;
# Defining particle names.
POSSIBLE_SECTORS = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
# Tracking df rows.
rows = [];

# Purpose: To analyze stats pertaining to contact pairs and
# particle properties, particle types, sector interactions,
# sequence separation, and number of shared contacts.
# Parameters:
#   edgeSet = The set of edges describing the structure.
#   particleContacts = A dictionary containing the contact set for
#                      each particle in the structure.
#   parIDMap = A dictionary mapping particles to their IDs.
#   structID = The ID of the current structure.
def analyze(edgeSet, particleContacts, parIDMap, structID):
    for edge in edgeSet:
        permanentStatus = 0;
        try:
            resDiff = abs(edge[i_EDGE_RES_A] - edge[i_EDGE_RES_B]);
            parA_isBB = edge[i_EDGE_NAME_A][-1] == "0";
            parB_isBB = edge[i_EDGE_NAME_B][-1] == "0";
            # Only happens with R group bonds or peptide bonds.
            if resDiff == 0 or (resDiff == 1 and parA_isBB and parB_isBB):
                permanentStatus = 1;
        except:
            pass;
        analyzeEdge(edge, particleContacts, parIDMap, structID, permanentStatus);

# ------ Methods to Populate the Data Frame ------
# Purpose: Populating with real and fake contact information.
# Parameters:
#   edge = The current edge being analyzed.
#   particleContacts = A dictionary of contact sets in the structure.
#   parIDMap = A dictionary mapping particles to their IDs.
#   structID = The ID of the current structure.
#   permanentStats = Indicating whether the interaction is permanent or not.
def analyzeEdge(edge, particleContacts, parIDMap, structID, permanentStatus):
    # Particle A information.
    A_type = edge[i_EDGE_NAME_A]; # Type - A0, G1, L0, etc.
    A_res = edge[i_EDGE_RES_A];
    A_sector = edge[i_EDGE_SECTOR_A];
    A_ss = edge[i_EDGE_SS_A];
    A_par = (A_type, A_res);
    # Particle B information.
    B_type = edge[i_EDGE_NAME_B];
    B_res = edge[i_EDGE_RES_B];
    B_sector = edge[i_EDGE_SECTOR_B];
    B_ss = edge[i_EDGE_SS_B];
    B_par = (B_type, B_res);
    # Building and adding row. Edge_Status will always be 1 as we are
    # only adding real edges here.
    rows.append(buildRow(A_par, B_par, A_sector, B_sector, A_ss, B_ss,
                         parIDMap, structID, permanentStatus, edgeStatus = 1));
    # Now building and adding a fake row.
    # Finding empty sectors which could be used for the
    # hard negative.
    if A_res == "-":
        A_emptySectors = {-1};
    else:
        A_contactSet = particleContacts[A_par];
        A_occupiedSectors = set([x[i_SECTOR_NUM] for x in A_contactSet]);
        A_emptySectors = POSSIBLE_SECTORS - A_occupiedSectors;
    if B_res == "-":
        B_emptySectors = {-1};
    else:
        B_contactSet = particleContacts[B_par];
        B_occupiedSectors = set([y[i_SECTOR_NUM] for y in B_contactSet]);
        B_emptySectors = POSSIBLE_SECTORS - B_occupiedSectors;
    # If there are no empty sectors, no hard negative can be made.
    if len(A_emptySectors) > 0 and len(B_emptySectors) > 0:
        A_newSector = rand.choice(list(A_emptySectors));
        B_newSector = rand.choice(list(B_emptySectors));
        rows.append(buildRow(A_par, B_par, A_newSector, B_newSector, A_ss, B_ss,
                             parIDMap, structID, permanentStatus, edgeStatus = 0));
    # Only making an easy negative for non-permanent entries.
    if not permanentStatus:
        while True: # Looping until we are able to make the easy negative
            # Picking two random particles.
            rand_A_par = rand.choice(list(particleContacts.keys()));
            rand_B_par = rand.choice(list(particleContacts.keys()));
            while rand_B_par == rand_A_par: # Making sure we don't pick the same particle twice.
                rand_B_par = rand.choice(list(particleContacts.keys()));
            # Picking random empty sectors for the edge..
            rand_A_contactSet = particleContacts[rand_A_par];
            rand_B_contactSet = particleContacts[rand_B_par];
            # Finding empty sectors for both particles. If they don't exist, we should
            # go back to selection.
            rand_A_occupiedSectors = set([x[i_SECTOR_NUM] for x in rand_A_contactSet]);
            rand_A_emptySectors = POSSIBLE_SECTORS - rand_A_occupiedSectors;
            rand_B_occupiedSectors = set([y[i_SECTOR_NUM] for y in rand_B_contactSet]);
            rand_B_emptySectors = POSSIBLE_SECTORS - rand_B_occupiedSectors;
            # If there are no empty sectors on either, we should redo everything.
            if len(rand_A_emptySectors) == 0 or len(rand_B_emptySectors) == 0:
                continue;
            # Picking random sectors.
            rand_A_sector = rand.choice(list(rand_A_emptySectors));
            rand_B_sector = rand.choice(list(rand_B_emptySectors));
            # Information pertaining to both particles.
            rows.append(buildRow(rand_A_par, rand_B_par, rand_A_sector, rand_B_sector, A_ss, B_ss,
                                 parIDMap, structID, permanentStatus, edgeStatus = 0));
            break; # Loop will end.

# Purpose: To build a df row dataframe.
# Parameters:
#   A_par = A particle.
#   B_par = B_particle.
#   A_sector = The sector for particle A.
#   B_sector = The sector for particle B.
#   A_ss = The secondary structure value for particle A.
#   B_ss = The secondary structure value for particle B.
#   parIDMap = A dictionary mapping all particles to their IDs.
#   structID = The structure ID.
#   permanentStatus = The permanent status of the edge.
#   edgeStatus = 1 or 0 indicating if the edge exists.
# Return:
#   A df representing a row of the df.
def buildRow(A_par, B_par, A_sector, B_sector, A_ss, B_ss, parIDMap, structID, permanentStatus, edgeStatus):
    # Par A information.
    A_type = A_par[0];  # Type - A0, G1, L0, etc.
    A_res = A_par[1];
    A_ID = parIDMap[A_par];
    # Par B information.
    B_type = B_par[0];  # Type - A0, G1, L0, etc.
    B_res = B_par[1];
    B_ID = parIDMap[B_par];
    # Information pertaining to both particles.
    try:
        seqDiff = A_res - B_res;
        seqSep = abs(seqDiff);
        solventStatus = 0;
    except:
        seqDiff = 0.5;
        seqSep = 0.5;
        solventStatus = 1;
    # Building and returning row.
    return (structID, A_ID, A_type, A_sector, A_ss, B_ID, B_type, B_sector, B_ss,
            seqSep, seqDiff, solventStatus, permanentStatus, edgeStatus);

# Purpose: To export the dataframe.
# Parameters:
#   outputDir = The output directory.
#   outputSuff = The suffix for the output file.
def export(outputDir, outputSuff):
    df = pd.DataFrame(rows, columns = ["Struct_ID",
                                       "ID_A", "Type_A", "Sector_A", "SS_A",
                                       "ID_B", "Type_B", "Sector_B", "SS_B",
                                       "Sequence_Separation", "Sequence_Difference",
                                       "Solvent_Status", "Permanent_Status", "Edge_Status"]);
    df.to_csv(f"{outputDir}/edges_{outputSuff}.tsv", sep = "\t", index = False);