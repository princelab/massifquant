require 'spec_helper'

SINGLE_LC_PEAK =<<END
794.153320312500	66.6232154134822
794.155700683594	64.1848374569493
794.152343750000	77.3453598317313
794.137878417969	43.9971354158960
794.153381347656	47.9325800115381
794.151916503906	58.2493545430977
794.151794433594	48.2357618443321
794.154418945313	73.7836456620808
794.155273437500	64.6420998592926
794.154724121094	50.6784271511619
794.153137207031	75.8204119171167
794.152648925781	46.7605242791329
794.154663085938	52.1323014595438
794.154479980469	84.6594395393981
794.153686523438	100.516711918964
794.153442382813	100.083305145014
794.153808593750	123.268559232728
794.154052734375	103.280867056779
794.152282714844	81.9609252126814
794.154174804688	81.7546736288284
794.154846191406	105.340162353622
794.154235839844	90.3293966676616
794.153381347656	112.033485856350
794.154296875000	118.121734412533
794.153442382813	110.789179018079
794.153076171875	122.038260041933
794.154724121094	73.9127937113647
794.154174804688	113.927951464950
794.153747558594	145.712131926703
794.154357910156	125.156195380762
794.152770996094	73.7180247155928
794.153869628906	125.161871752643
794.153564453125	121.757038265504
794.153930664063	131.414978158409
794.154418945313	172.332238396143
794.153808593750	150.444848179491
794.153442382813	142.661866674753
794.154174804688	155.523172007421
794.153442382813	142.257661261178
794.153564453125	146.678159321267
794.154418945313	146.469673428922
794.154296875000	156.684340065863
794.153503417969	169.645966371883
794.153625488281	174.528007467641
794.153259277344	131.354965468097
794.153686523438	130.782343932103
794.154357910156	131.504076992407
794.153686523438	156.862985097664
794.153869628906	113.358925955348
794.155517578125	64.6764219076879
794.155822753906	77.2962606732208
794.155639648438	82.3411793953594
794.154296875000	76.4721852518123
794.154785156250	80.5094132695053
794.152526855469	54.2418426146677
794.153747558594	49.0973328776943
END


# ghost scan properties
[:IntVar, :MZVar, :MaxIntSqrt, :MaxScanNum, :MaxDataNum, :I, :MZ]
[4.936164570985432e+07, 5.774199962615967e-07, 1.766446600877026e+02, 18, 4, [1.268376464843750e+04, 2.147245703125000e+04, 3.120333593750000e+04, 2.419581250000000e+04, 1.720981640625000e+04], [7.944052734375000e+02, 7.944046020507812e+02, 7.944047851562500e+02, 7.944033203125000e+02, 7.944050292968750e+02]]

describe "Kalmanquant on a single chromatographic peak" do
  before do
    # this simulates the data structure of an entire run (hence the nesting)
    @ar_of_spectrum = SINGLE_LC_PEAK.each_line.map {|line| line.split(/\s+/).map {|v| [v.to_f] } }
  end

  xit "tracks a single peak" do
    peaks = Kalmanquant.find_features(@ar_of_spectrum)
    peaks.nil?.is false
    p @mzs 
    p @intensities
    ##1.is 1
  end
end
