############################################################################
############################################################################
# Jae Wook Cho                                                             #
# 08/12/19 Copy-Paste translate from:                                      #
# https://github.com/ebridge2/Discriminability/blob/master/R/reliability.R #
############################################################################
############################################################################
# Reliability Density function:
def discr_rdf(dist, ids):
    """
    Outputs the reliability density function of the entire subject similarity matrix.

    Dist is a subject similarity matrix.
    ids is subject ids as strings.
    """
    import numpy as np
    import sys
    N = dist.shape[0]
    if N <= 1:
        sys.exit('Invalid datatype for N')

    uniqueids = np.unique(ids)
    countvec = np.zeros(int(len(uniqueids)))

    for i in range(len(uniqueids)):
        countvec[i] = np.sum(uniqueids[i] == ids) # total number of scans for particular id

    scans = np.max(countvec) # assume that we will worst case have the most
    rdf = np.zeros(int(N*(scans-1))) # initialize empty ra

    count = 0
    for i in range(N):
        ind = np.where(ids[i] == ids)[0] # all the indices that are the same subject, but different scan
        for j in ind:
            if i != j:
                di = dist[i,:].copy()
                di[ind] = float('Inf')
                d = dist[i,j]
                rdf[count] = 1 - (np.nansum(di < d) + (0.5*np.nansum(di == d))) / (N - len(ind))
                count = count + 1
    return rdf[0:count]


# Cross condition Reliability Density function:
def discr_xrdf(dist, sourceinds, targetinds,targetids):
    """
    Outputs the reliability density function between 2 conditions.

    Dist is a subject similarity matrix.
    sourceinds are the indices of the source condition.
    targetinds are the indices of the target condition.
    targetids are subject ids as strings.
    """
    import numpy as np
    import sys
    dist = dist[sourceinds,:][:,targetinds]
    ids = targetids
    N = dist.shape[0]
    if N <= 1:
        sys.exit('Invalid datatype for N')

    uniqueids = np.unique(ids)
    countvec = np.zeros(int(len(uniqueids)))

    for i in range(len(uniqueids)):
        countvec[i] = np.sum(uniqueids[i] == ids) # total number of scans for particular id

    scans = np.max(countvec) # assume that we will worst case have the most
    rdf = np.zeros(int(N*(scans-1))) # initialize empty ra

    count = 0
    for i in range(N):
        ind = np.where(ids[i] == ids)[0] # all the indices that are the same subject, but different scan
        for j in ind:
            if i != j:
                di = dist[i,:].copy()
                di[ind] = float('Inf')
                d = dist[i,j]
                rdf[count] = 1 - (np.nansum(di < d) + (0.5*np.nansum(di == d))) / (N - len(ind))
                count = count + 1
    return rdf[0:count]

# Discriminability
def discr_discr(rdf, remove_outliers=True, thresh=0, output=False):
    """
    Discriminability function:
    
    Input is the output from either discr_rdf or discr_xrdf.
    Outputs 1 discriminability value.
    """
    import numpy as np
    if remove_outliers:
        discr = np.mean(rdf[np.where(rdf > thresh)[0]]) # mean of the rdf
        ol = len(np.where(rdf < thresh)[0])
        if output:
            print('Graphs with reliability < %s (outliers): %s' % (thresh, ol))

    else:
        ol = 0
        discr = np.nanmean(rdf)

    nopair = np.sum(np.isnan(rdf))
    if output:
        print('Graphs with unique ids: %s' % nopair)
        print('Graphs available for reliability analysis: %s' % (len(rdf)-ol-nopair))
        print('discr: %s' % discr)

    return(discr)