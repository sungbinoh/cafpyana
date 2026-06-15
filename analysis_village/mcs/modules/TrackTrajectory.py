import numpy as np

class Trajectory():
    """
    A class to represent a track trajectory.
    """
    def __init__(self, df):
        """
        Initialize the trajectory with a dataframe of hit information.
        Dataframe must have these columns: x, y, z, rr, pitch, tpc
        """
        self.hits = df.copy() 

        # TODO: cut on pitch?
        valid = (self.hits.rr > 0)
        self.hits = self.hits[valid]

        # sort by rr in descending order
        self.hits.sort_values("rr", inplace=True, ascending=False)

        self.nhits = len(self.hits)

        # trajectory length is the sum of Euclidean distances between consecutive hits
        distances = np.linalg.norm(self.hits[["x", "y", "z"]].diff(), axis=1)[1:] # first entry will be nan
        #distances = np.linalg.norm(self.hits[[('h', 'sp', 'x', ''), ('h', 'sp', 'y', ''), ('h', 'sp', 'z', '')]].diff(), axis=1)[1:] # first entry will be nan
        
        self.length = np.sum(distances)

    def first_point(self):
        """
        Return the index of the first hit in the trajectory.
        """
        return self.hits.index[0]
    
    def next_point(self, idx):
        """
        Return the index of the next hit in the trajectory.
        """
        hits = self.hits
        if idx == hits.index[-1]:
            return None # last hit
        return hits.index[hits.index.get_loc(idx)+1]

    def last_point(self):
        """
        Return the index of the last hit in the trajectory.
        """
        return self.hits.index[-1]
    

    def location_at_point(self, idx):
        """
        Return the location of the hit at the given index.
        """
        x = self.hits.loc[idx].x
        y = self.hits.loc[idx].y
        z = self.hits.loc[idx].z
        
        
        #print("[location_at_point] x")
        #print(x)
        #print(type(x))
        return np.array([float(x),float(y),float(z)])

    
    def direction_at_point(self, idx):
        """
        Return the direction of the trajectory at the given index.
        """
        next_idx = self.next_point(idx)
        pos_this = self.location_at_point(idx)
        pos_next = self.location_at_point(next_idx)
        #print("[direction_at_point] pos_this")
        #print(pos_this)
        vec = pos_next-pos_this
        mag = np.linalg.norm(vec)
        if mag == 0:
            return None
        norm_vec = vec/mag
        return norm_vec


    def rr_at_point(self, idx):
        """
        Return the residual range of the hit at the given index.
        """
        return float(self.hits.loc[idx].rr)


    def pitch_at_point(self, idx):
        """
        Return the pitch of the hit at the given index.
        """
        return float(self.hits.loc[idx].pitch)


    def tpc_at_point(self, idx):
        """
        Return the TPC of the hit exists in at the given index.
        """
        return self.hits.loc[idx].tpc
