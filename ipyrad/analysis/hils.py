


class Hils(object):
    """
    A Class to calculate the Hils statistic given a matrix of invariants.
    """
    def __init__(self, database, boot=0, tree=None, root=None):
        ## open file handles for accessing database
        self._open = True
        self._boot = boot
        self.hdf5 = h5py.File(database)
        self.matrix = self.hdf5["invariants"]["boot{}".format(self._boot)]
        self.quartets = self.hdf5["quartets"]
        self.nquartets = self.quartets.shape[0]
        self.tree = tree
        self.root = root
        if self.tree:
            self.snames = sorted(self.tree.get_tip_labels())
            self.sidx = {i:j for i,j in enumerate(snames)}
        

    def close_db(self):
        """close the database file"""
        self.hdf5.close()
    
    

    def get_counts_by_idx(self, idx, altmat=None):
        """
        Return site counts for a given index (quartet). Chooses the 
        'correct' matrix based on the name order in self.quartets. 
        But this can be overridden during testing by entering a 
        altmat index.
        """

        ## the matrix is stored in default order format (e.g., 0,1|2,3)
        mat = self.matrix[idx, :, :]

        ## the correct quartet is stored separate (e.g., 0,3|1,2)
        qrt = self.quartets[idx]
        
        ## the matrix needs to be arranged to be in the right order.
        ## if taxon 1 is the second lowest (e.g., 0,1|2,3) then no reorder
        ## if taxon 1 is the third lowest (e.g., 0,2|1,3) then reorder mat1
        ## if taxon 1 is the highest (e.g., 0,3|1,2) then reorder to mat2
        if isinstance(altmat, int):
            assert altmat in [0, 1, 2], "altmat must be an index in [0,1,2]"
            mat = alt_mats(mat, altmat)
        else:
            if qrt[1] > qrt[2]:
                if qrt[1] > qrt[3]:
                    mat = alt_mats(mat, 2)
                else:
                    mat = alt_mats(mat, 1)
            
        ## return counts as a dataframe with column names
        df = pd.DataFrame(
            data=count_snps(mat), 
            index=["aabb", "abba", "baba", "aaab"], 
            columns=[idx]).T
        return df
    

    
    def get_h_by_idx(self, idx, altmat=None):
        """
        calculate Hils. This could be numba-fied, but you'd have to work
        with arrays instead of dataframes. This is fine for now.
        """

        ## get counts and convert to site frequencies
        df = self.get_counts_by_idx(idx, altmat)
        nsites = df.sum(axis=1).values[0]
        pdf = df/nsites
        pdf.columns = ["p"+i for i in df.columns]
        data = pd.concat([df, pdf], axis=1)

        ## avoid zero div errors
        if data.pabba.equals(data.pbaba):
            H = 0.0
            f1 = 1.0
            f2 = 0.0

        else:
            ## get H and f1 and f2 for these data
            H, f1, f2 = calc_h(data, nsites)

            ## f1 and f2 measure differences/distances, should be positive
            f1, f2 = [abs(i) for i in (f1, f2)]

        ## return as a dataframe 
        res = pd.DataFrame(
             {"Hils":H,
              "gamma": 1. - (f1/(f1+f2)),
              "pval": norm.pdf(H, 0, 1)}, 
             index=[idx],
             )
        return pd.concat([df, pdf, res], axis=1)



    def run(self):
        """calculate Hils and return table for all idxs in database"""
        stats = pd.concat([self.get_h_by_idx(idx) for idx in xrange(self.nquartets)])
        qrts = ["{},{}|{},{}".format(*i) for i in self.quartets[:]]
        qrts = pd.DataFrame(np.array(qrts), columns=["qrts"])
        return pd.concat([stats, qrts], axis=1)




    def svds(self, idx):
        """
        returns the svd scores for the three resolutions of the matrix
        as calculated by tetrad. 
        """
        mats = np.zeros((3, 16, 16), dtype=np.uint32)
        mats[0] = self.matrix[idx]
        mats[1] = alt_mats(mats[0], 1)
        mats[2] = alt_mats(mats[0], 2)

        svds = np.zeros((3, 16), dtype=np.float64)
        scor = np.zeros(3, dtype=np.float64)
        rank = np.zeros(3, dtype=np.float64)

        ## why svd and rank?
        for test in range(3):
            svds[test] = np.linalg.svd(mats[test].astype(np.float64))[1]
            rank[test] = np.linalg.matrix_rank(mats[test].astype(np.float64))

        ## get minrank, or 11
        minrank = int(min(11, rank.min()))
        for test in range(3):
            scor[test] = np.sqrt(np.sum(svds[test, minrank:]**2))

        ## sort to find the best qorder
        return scor

    

def calc_h(data, nsites):
    """ 
    Calculate Hils statistic from site counts/frequencies.
    """

    f1 = data.paabb - data.pbaba
    f2 = data.pabba - data.pbaba           

    sigmaf1 = (1. / nsites) * (data.paabb * (1. - data.paabb) \
        + data.pbaba * (1. - data.pbaba) \
        + 2. * data.paabb * data.pbaba)

    sigmaf2 = (1. / nsites) * (data.pabba * (1. - data.pabba) \
        + data.pbaba * (1. - data.pbaba) \
        + 2. * data.pabba * data.pbaba)

    covf1f2 = (1. / nsites) * (data.pabba * (1. - data.paabb) \
        + data.paabb * data.pbaba \
        + data.pabba * data.pbaba \
        + data.pbaba * (1. - data.pbaba)) 

    num = f2 * ((f1 / f2) - 0.)
    p1 = (sigmaf2 * (f1/f2)**2)
    p2 = ((2. * covf1f2 * (f1/f2) + sigmaf1))
    denom = p1 - p2

    ## calculate hils
    H = num/np.sqrt(abs(denom))
    return H, f1, f2

    

@numba.jit(nopython=True)   
def alt_mats(mat, idx):
    """ return alternate rearrangements of matrix"""
    mats = np.zeros((3, 16, 16), dtype=np.uint32)
    mats[0] = mat
    x = np.uint8(0)
    for y in np.array([0, 4, 8, 12], dtype=np.uint8):
        for z in np.array([0, 4, 8, 12], dtype=np.uint8):
            mats[1, y:y+np.uint8(4), z:z+np.uint8(4)] = mats[0, x].reshape(4, 4)
            mats[2, y:y+np.uint8(4), z:z+np.uint8(4)] = mats[0, x].reshape(4, 4).T
            x += np.uint8(1)
    return mats[idx]
        
        

@numba.jit(nopython=True)
def count_snps(mat):
    """ JIT func to return counts quickly."""
    ## array to store results
    snps = np.zeros(4, dtype=np.uint16)

    ## get concordant (aabb) pis sites
    snps[0] = np.uint16(\
           mat[0, 5] + mat[0, 10] + mat[0, 15] + \
           mat[5, 0] + mat[5, 10] + mat[5, 15] + \
           mat[10, 0] + mat[10, 5] + mat[10, 15] + \
           mat[15, 0] + mat[15, 5] + mat[15, 10])

    ## get discordant (baba) sites
    for i in range(16):
        if i % 5:
            snps[1] += mat[i, i]

    ## get discordant (abba) sites
    snps[2] = mat[1, 4] + mat[2, 8] + mat[3, 12] +\
              mat[4, 1] + mat[6, 9] + mat[7, 13] +\
              mat[8, 2] + mat[9, 6] + mat[11, 14] +\
              mat[12, 3] + mat[13, 7] + mat[14, 11]

    ## get autapomorphy sites
    snps[3] = (mat.sum() - (snps[0] + np.diag(mat).sum() + snps[2]))
    return snps

