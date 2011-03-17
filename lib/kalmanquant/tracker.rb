require 'matrix'

class Matrix
  # A\B == A^-1 * B
  def backslash(other)
    self.inverse * other
  end

  def []=(i, j, x)
    @rows[i][j] = x
  end
end

module Kalmanquant
  # finds a large peak and sets reasonable estimates for variances
  GhostScan = Struct.new(:mz_var, :int_var, :max_int) do
    # returns a GhostScan object based on the spectra
    def self.new(spectra)
    end

    def self.default_vals
      self.new(5.7742e-07, 4.9362e07, 31203.35003809)
    end
  end
  class Tracker

    MZQ = Matrix.zero(2)
    G = Matrix[[0, 0]]
    U = MZG.transpose
    G_X_U = G * U

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

    # the total number of spectra seen (sum of updates and misses)
    attr_reader :num_spectra_total

    # a ghost scan Chris found:
    # m/zs (just 5 values):   [794.4053 794.4046 794.4048 794.4033 794.4050]
    # intensities:  [   1.2684 2.1472 3.1203 2.4196 1.7210 ] * 10e4
    # IntVar: 4.9362e+07
    # MZVar: 5.7742e-07  # could square root this (which would make it larger)
    # MaxIntSqrt: 176.6447

    def initialize(peak, current_time, ghost_scan=GhostScan.default_vals)
      @times << current_time
      @peaks = [peak] 
      @mzr = Math.sqrt(ghost_scan.mz_var)
      @intq = Matrix.identity(2) * ghost_scan.int_var
      @intr = Math.sqrt(ghost_scan.max_int)
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

    # returns self
    def prediction!(time)
      # this is not in Chris's code
      delta_time = (time - @times.last).abs
      @mzf[0,1] = delta_time
      @intf[0,1] = delta_time
      # prediction mass
      @mzp = @mzf*@mzp*@mzf.transpose + MZQ
      # the bit on the end about mzg and mzu is not really doing anything ......
      # tr.MZXhat = tr.MZF*tr.MZXhat + tr.MZG*tr.MZU;
      @mz_pred = @mzf*@mz_pred + G_X_U

      # prediction intensity
      #tr.IntP = tr.IntF*tr.IntP*tr.IntF' + tr.IntQ;
      @intp  = @intf*@intp*@intf.transpose + @intq
      # the bit on the end is not doing anything
      # tr.IntXhat = tr.IntF*tr.IntXhat + tr.IntG*tr.IntU;
      @int_pred = @intf*@int_pred+ G_X_U
      self
    end

    def mz_margin ; Math.sqrt(3*@mzp[0,0]) end
    def int_margin ; Math.sqrt(3*@intp[0,0]) end

    # returns the lo and hi predictions
    def pred_mz_margins
      mz_m = mz_margin
      mz_pred_00 = @mz_pre[0,0]
      [mz_pred_00-mz_m, mz_pred_00+mz_m]
    end

    def pred_int_margins
      int_m = int_margin
      int_pred_00 = @int_pre[0,0]
      [int_pred_00-int_m, int_pred_00+int_m]
    end

    def missed! 
      @num_spectra_total += 1
      @num_missed += 1 
    end

    def normalized_distance(mz, int)
      # Cinf=(m/z(p)-m/z(m))^2/std(m/z)+(I(p)-I(m))^2/std(I);
      ((mz - @mz_pred)**2 / Math.sqrt(@mzp[0,0])) + ((int - @int_pred)**2 / Math.sqrt(@intp[0,0]))
    end

    # this should be deprecated in favor of asking if normalized_distance is
    # close enough.
    def in_bounds?(mz, int)
      (pred_lo_mz, pred_hi_mz) = pred_mz_margins
      (pred_lo_int, pred_hi_int) = pred_int_margins
      (mz >= pred_lo_mz) && (mz <= pred_hi_mz) && (int >= pred_lo_int) && (int <= pred_hi_int)
    end

    # updates current position (i.e., the innovation step) based on the new peak
    # [mzp, mz_pred, intp, int_pred].  Sets num_missed back to zero. returns
    # self
    def update!(peak, current_time)
      @num_spectra_total += 1
      # Ralf and Chris do this slightly differently conceptually (although it
      # should be equivalent based on Jeff's paper)
      @peaks << peak
      @times << current_time
      @num_missed = 0
      @mzp = (@mzp.inverse + H_TRANSFORM_TRANSPOSE * @mzr.backslash(H_TRANSFORM)).inverse
      @mz_pred = @mz_pred - @mzp*(H_TRANSFORM_TRANSPOSE/@mzr) * (H_TRANSFORM*@mz_pred - peak.mz)
      @int_pred = @int_pred - @intp*(H_TRANSFORM_TRANSPOSE/@intr) * (H_TRANSFORM*@int_pred - peak.intensity)
      self
    end

  end
end

