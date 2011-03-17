require 'spec_helper'

require 'kalmanquant/tracker'
require 'kalmanquant/tracker_manager'

FakeTracker = Struct.new(:lo_mz, :hi_mz, :pred_mz, :lo_int, :hi_int, :pred_int, :size, :miss_cnt) do
  def distance(mz, int)
    (pred_mz - mz).abs
  end
  def in_bounds?(mz, int)
    mz >= self.lo_mz && mz <= self.hi_mz
  end
  def missed!
    @miss_cnt += 1
  end
  def update!(*args)
    @update_data = args
  end
end

FakeSpectrum = Struct.new(:mzs, :intensities)

describe 'finding matches between spectra and trackers' do
  before do
    @spectrum = FakeSpectrum.new([1.1, 1.4, 1.8, 1.9, 2.2, 3.6, 4.5], Array.new(7) {500} )
    # 
    # |-------------X-|       # 0 -> 0
    # |--------------X-----|  # 1 -> 1
    # |---------------------------------X-----------------------------------------| # 2 -> 3
    #                                              |---X---|  # 3 -> nil
    #                                                |---X---|  # 4 -> 5
    #                                                              |---X---| # 5-> nil
    #                                                                    |---X---|  # 6 -> 6 (this guy has large size, so wins!)
    #           1.1    1.4      1.8 1.9    2.2             3.6           4.5
    @trackers = [[0,1.5, 1.2], [0,1.5,1.21], [0,5,1.95], 
      [2.6, 3.5, 3.2], [2.7, 3.6, 3.4],
      [4.2, 4.8, 4.4], [4.4, 4.8, 4.6]].map {|row| FakeTracker.new(*row, 0, 1000, 500, 1) }
    @tm = Kalmanquant::TrackerManager.new
    @best_matches_key = [[0,0], [1,1], [2,3], [4,5], [6,6]]
  end

  it 'works on very sparse data' do
    @tm.possible_matches(FakeSpectrum.new([], []), @trackers).empty?.is true
    @tm.possible_matches(@spectrum, []).empty?.is true
    @tm.possible_matches(FakeSpectrum.new([1.1], [1000]), [FakeTracker.new(1.3, 1.5, 1.4, 999, 1001, 1000)]).empty?.is true
  end

  it 'finds matches greedily' do
    @trackers[-1].size = 2
    (best_matches, tracker_set, index_set) = @tm.best_matches(@spectrum, @trackers)
    best_matches.size.is @best_matches_key.size
    @best_matches_key.each do |track_i, mz_i|
      matches = best_matches.select {|match| match.tracker == @trackers[track_i] }
      matches.size.is 1
      matches.first.centroid_index.is mz_i
    end

    # now, make the last tracker size smaller than the others
    @trackers[-1].size = 0
    (another_try, _, _) = @tm.best_matches(@spectrum, @trackers)
    another_try.select {|match| match.tracker == @trackers[-1] }.empty?.is true
    another_try.select {|match| match.tracker == @trackers[-2] }.size.is 1

    @trackers[-1].size = 1
  end

  it 'finds possible matches' do
    matches = @tm.possible_matches(@spectrum, @trackers)
  
    match0 = matches.select {|match| match.tracker == @trackers[0] }
    match0.size.is 2

    match2 = matches.select {|match| match.tracker == @trackers[2] }
    match2.size.is 7
    # frozen, but checked for basic sanity:
    match2.sort_by(&:centroid_index).map(&:distance).zip([0.85, 0.55, 0.15, 0.05, 0.25, 1.65, 2.55]) do |act, exp| 
      act.should.be.close exp, 0.00000001 
    end
    matches.select {|match| match.tracker == @trackers[3] }.empty?.is true
    match4 = matches.select {|match| match.tracker == @trackers[4] }
    match4.size.is 1
    match4.first.distance.should.be.close 0.2, 0.000000001
    match4.first.centroid_index.is 5
  
    matches.empty?.is false
  end
end
