from Core.fetchInfo import i_PAR_NAME, i_PAR_AA, i_PAR_TYPE, i_PAR_RES, i_PAR_SS,\
    i_CONTACT_NAME, i_CONTACT_AA, i_CONTACT_TYPE, i_CONTACT_RES, i_CONTACT_SS,\
    i_CONTACT_AREA, i_SECTOR_NUM, i_RES_DIFF;

# Index constants for keys.
i_KEY_NAME = 0;
i_KEY_RES = 1;

# Purpose: Populating the contents of a given file's
# contact sets into an edge set and returning it.
# Parameters:
#   contactSets = A dictionary with keys representing
#                 particles and values representing
#                 their contact sets as arrays.
# Return:
#   edgeSet = The set of edges in the given contact set.
def getEdgeSet(contactSets):
    edgeSet = set();
    # Iterating through the particles.
    for parKey in contactSets.keys():
        parContacts = contactSets[parKey];
        parName = parKey[i_KEY_NAME];
        parRes = parKey[i_KEY_RES];
        # Iterating through the contacts, making the edges.
        for contact in parContacts:
            parSector = contact[i_SECTOR_NUM];
            parSS = contact[i_PAR_SS];
            # Collecting information about the contact.
            contactName = contact[i_CONTACT_NAME];
            contactRes = contact[i_CONTACT_RES];
            contactSS = contact[i_CONTACT_SS];
            if contactName == "O0": # If a solvent molecule, assume the corresponding edge exists.
                edge = (parName, parRes, parSector, parSS, contactName, contactRes, -1, contactSS); # Solvent sector is always -1.
                edgeSet.add(edge); # Cannot have reverse edges.
            else:
                contactKey = (contactName, contactRes);
                # Now that we have the contact key, we can
                # find the sector corresponding to this
                # contact.
                corrSet = contactSets[contactKey];
                corrContact = [x for x in corrSet if x[i_CONTACT_NAME] == parName and x[i_CONTACT_RES] == parRes];
                try:
                    contactSector = corrContact[0][i_SECTOR_NUM];
                    # Now that we have the corresponding sector,
                    # we can make the edge.
                    edge = (parName, parRes, parSector, parSS, contactName, contactRes, contactSector, contactSS);
                    reverseEdge = (contactName, contactRes, contactSector, contactSS, parName, parRes, parSector, parSS);
                    # Making sure not to add redundant information.
                    if reverseEdge not in edgeSet:
                        edgeSet.add(edge);
                # No corresponding contact exists, so this edge does
                # not actually exist. Just pass, removing the contact
                # from the current contact set.
                except:
                    parContacts.remove(contact);
    return edgeSet;