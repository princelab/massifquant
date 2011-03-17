class Matrix
  # A\B == A^-1 * B
  def backslash(other)
    self.inverse * other
  end
end


module Kalmanquant
  class Tracker

    MZQ = Matrix.identity(2) * 10e-4
    # throws away velocity component (a column vector) (MZH in chris's code)
    H_TRANSFORM = Matrix[[1, 0]]
    H_TRANSFORM_TRANSPOSE = H_TRANSFORM.transpose

    # the array of spectrum indices that this tracker tracks
    attr_accessor :peaks

    attr_accessor :max_int_sqrt

    # Integer, initialized to zero
    attr_accessor :num_missed

    # a 2x1 Matrix
    attr_accessor :mz_pred
    # a 2x1 Matrix
    attr_accessor :int_pred

    # a 2X2 matrix
    attr_accessor :mzp
    # a 2X2 matrix
    attr_accessor :intp

    # measurement uncertainty (can be a matrix in other implementations; here
    # it is a scalar)
    attr_accessor :mzr

    # intensity uncertainty (can be a matrix in other implementations; here
    # it is a scalar)
    attr_accessor :intr

    # the int process variance (a 2X2 matrix)
    attr_accessor :intq
     
    # tr.MZR = ghost.MZVar
    # tr.IntQ = ghost.IntVar.*tr.IntQ;
    # tr.IntR = ghost.MaxIntSqrt*tr.IntR;

    # returns MZQ (the m/z process variance (a 2X2 matrix))
    def mzq
      MZQ
    end

    # a ghost scan Chris found:
    # m/zs (just 5 values):   [794.4053 794.4046 794.4048 794.4033 794.4050]
    # intensities:  [   1.2684 2.1472 3.1203 2.4196 1.7210 ] * 10e4
    # IntVar: 4.9362e+07
    # MZVar: 5.7742e-07  # could square root this (which would make it larger)
    # MaxIntSqrt: 176.6447


    def initialize(peak, delta_time, ghost_mz_var=0.1, ghost_int_var=1000, ghost_max_int_sqrt=1)
      @peaks = [peak] 
      @mzr = ghost_mz_var
      @intq = Matrix.identity(2) * ghost_int_var
      @intr = ghost_max_int_sqrt
      @mz_pred = Matrix[[peak.mz],[0]]
      @int_pred = peak.intensity
      @spectrum_ids = [spectrum_id]
      @centroid_ids = [centroid_id]
      # use mz_var to set mzp
      @mzp = Matrix.identity(2)
      @mzp[0,0] = mz_var
      # use int_var to set intp
      @intp = Matrix.identity(2)
      @intp[0,0] = int_var
      @max_int_sqrt = max_int_sqrt
    end

    # need to call prediction step on ALL active trackers  (NEED TO KEEP TRACK
    # of the TIMES so that the DT is accurate!!!!!!!!!!!!! (for trackers that
    # missed)] 
    def prediction!(delta_time)
      # precisely the same conceptually between Chris and Ralf
        # %prediction mass
        # tr.MZP = tr.MZF*tr.MZP*tr.MZF' + tr.MZQ;
        # tr.MZXhat = tr.MZF*tr.MZXhat + tr.MZG*tr.MZU;

        # %prediction intensity
        # tr.IntP = tr.IntF*tr.IntP*tr.IntF' + tr.IntQ;
        # tr.IntXhat = tr.IntF*tr.IntXhat + tr.IntG*tr.IntU;
    end

    def mz_margin ; Math.sqrt(3*@mzp[0,0]) end
    def int_margin ; Math.sqrt(3*@intp[0,0]) end

    # the low m/z cutoff for predicted peaks
    def pred_lo_mz ; @mz_pred[0,0] - mz_margin end
    # the high m/z cutoff for predicted peaks
    def pred_hi_mz ; @mz_pred[0,0] + mz_margin end
    def pred_lo_int ; @int_pred[0,0] - int_margin end
    def pred_hi_int ; @int_pred[0,0] + int_margin end

    def missed! ; @num_missed += 1 end

    # Cinf=(m/z(p)-m/z(m))^2/std(m/z)+(I(p)-I(m))^2/std(I);
    def normalized_distance(mz, int)
      ((mz - @mz_pred)**2 / Math.sqrt(@mzp[0,0])) + ((int - @int_pred)**2 / Math.sqrt(@intp[0,0]))
    end

    def in_bounds?(mz, int)
      (mz >= pred_lo_mz) && (mz <= pred_hi_mz) && (int >= pred_lo_int) && (int <= pred_hi_int)
    end

    # updates current position (i.e., the innovation step) based on the new peak
    # mzp, mz_pred, intp, int_pred
    def update!(peak, current_time)
      # Ralf and Chris do this slightly differently conceptually (although it
      # should be equivalent based on Jeff's paper)
      @peaks << peak
      @mzp = (@mzp.inverse + H_TRANSFORM_TRANSPOSE * @mzr.backslash(H_TRANSFORM)).inverse
      @mz_pred = @mz_pred - @mzp*(H_TRANSFORM_TRANSPOSE/@mzr) * (H_TRANSFORM*@mz_pred - peak.mz)
      @int_pred = @int_pred - @intp*(H_TRANSFORM_TRANSPOSE/@intr) * (H_TRANSFORM*@int_pred - peak.intensity)
    end

  end
end

