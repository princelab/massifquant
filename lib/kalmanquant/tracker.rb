
module Kalmanquant
  class Tracker
    include Enumerable
    attr_accessor :times
    attr_accessor :mzs
    attr_accessor :intensities

    def size ; intensities.size end
    alias_method :length, :size

    def each(&block)
      times.zip(mzs, intensities, &block)
    end

  end
end

=begin

classdef Tracker7 < handle
   %% Key Properties of Tracker
   % it should be a handle class because it will be passed around a lot. No
   % need to make more copies than necessary as it is modified by a
   % function. Each tracker instance should be unique. It also allows 
   % events & listener functionality, a plus for communication with other 
   % trects.

   %%
    % The following properties can be initialized in the constructor
    % Note, this initialization step will be fast for trackers
    % because the values are only called once upon first initilization
    % of one instance of Tracker Class (see mathworks: OO Properties)
   properties (SetAccess = public)
      %% Identity & Status
      ID
      DataPointNumbers = zeros(100,1);   %Each Scan has a corresponding DataPointNums
      %tracker path %never will reach 100
      ScanNumbers = zeros(100, 1);        %Tracker evolves over multiple time scans
      Status = TrackerStatus.ACTIVE; 
      %variables determining tracker status 
      TrackerLength = 1;    %all trackers start out with one data point
      CurrentMissed = 0; 
      
      %% DEBUG MODE
      %estimate variance
      TempPMZ = zeros(100,1);
      TempPInt = zeros(100,1);
      
      InnoPMZ = zeros(100,1);
      InnoPInt = zeros(100,1);
      
      %conf. intervals
      MZConf = zeros(100,2);
      IntConf = zeros(100,2);
      
      %sequence of hit - miss
      HitMiss = ones(1,1);
      
      %predictions
      IntPred = zeros(100,1);
      MZPred = zeros(100,1);
      %Kalman corrected estimates
      IntInno = zeros(100,1);
      MZInno = zeros(100,1);
      %% MODEL 
      %%Create the model within each tracker because each state-space will
      %have its own final converged representation.
      
      MZXhat = zeros(2,1);
      IntXhat = zeros(2,1);
      
      %Intensity Model
      IntF = [1 1
              0 1];
      IntG = [0 0];
      IntU = [0 
              0];
      IntH = [1 0]; 
      IntQ = eye(2);
      IntR = eye(1);
      %Arbitrary Initialization
      IntP = 10^(-3).*eye(2);
      
      %M/Z Model
      MZF = [1 1
             0 1];
      MZG = [0 0];
      MZU = [0 
             0];
      MZH = [1 0]; 
      MZQ = 10^(-4).*eye(2);
      MZR = eye(1);
      %Partially-Arbitrary Initialization
      MZP = 10^(1).*eye(2);
   end
   %%
   methods 
       %% Constructor
      function tr = Tracker7()
      end
      function copyCommonInit(tr, Init)
          %@ params: 
          %          Init -> [MZVar, IntVar,MaxIntSqrt, Scan #, dt];
          %
          % returns: An m x n array of Tracker trects.
          %necessary check for an array of trects
          if nargin > 0
             %key parameters are data dependent
             tr.MZR = Init(1);
             tr.IntQ = Init(2).*tr.IntQ;
             tr.IntR = Init(3)*tr.IntR;
             tr.ScanNumbers(1) = Init(4);
             %avg diff in time scans
             tr.MZF(1,2) = Init(5);
             tr.IntF(1,2) = Init(5);
             %initial variance
             tr.MZP(1,1) = Init(1);
             tr.IntP(1,1) = Init(2);
          end
      end
      function incrementMissed(tr)
          tr.CurrentMissed = tr.CurrentMissed + 1;
      end
      function resetMissed(tr)
          tr.CurrentMissed = 0;
      end
      function incrementTrackerLength(tr)
               tr.TrackerLength = tr.TrackerLength + 1;
      end
      function assignStatus(tr)
           %check tracker status
           if tr.CurrentMissed >= 2
              if tr.TrackerLength <= 2
                  tr.Status = TrackerStatus.NOISE;
              else
                  %tracker has claimed three or more centroids (real peak)
                  tr.Status = TrackerStatus.PIC;
              end
           else
                tr.Status = TrackerStatus.ACTIVE;
           end
      end
      function delete(h)
          %The DELETE method deletes a handle trect but does not clear the handle
          %from the workspace.  A deleted handle is no longer valid.
      end
      function predictCentroid(tr)
         %@ params: 
          %         tr -> The Tracker object
          %
          % returns: nothing, all info kept with the current trect
                     
         %prediction mass
         tr.MZP = tr.MZF*tr.MZP*tr.MZF' + tr.MZQ;
         tr.MZXhat = tr.MZF*tr.MZXhat + tr.MZG*tr.MZU;

         %prediction intensity
         tr.IntP = tr.IntF*tr.IntP*tr.IntF' + tr.IntQ;
         tr.IntXhat = tr.IntF*tr.IntXhat + tr.IntG*tr.IntU;

        
      end
      function [BestPointIdx LeastError] = claimCurrentCentroid(tr, CurrentScan) 
          %@ params: 
          %         tr -> The Tracker trect
          %          CurrentScan -> [m X 1] vector of data points
          % returns: best point claimed by tracker and its least error


           %                 ^    
           %i          
           %n          <   @ *     >
           %t              
           %
           %                 ^
           % % % % % % % % % % % % % % % % % %             
           %                  m/z
           
           %Marginal Error
           %evalsMZ = eigs(tr.MZP)
           %evalsInt = eigs(tr.IntP)
           MZMarginError = sqrt(3*tr.MZP(1,1));
           IntMarginError = sqrt(3*tr.IntP(1,1));
      
           %MZMarginError = 0.025;
           %IntMarginError = 10000;

           %Bounds on the Confidence Intervals
           Left = tr.MZXhat(1) - MZMarginError;
           Right = tr.MZXhat(1) + MZMarginError;
           Up =  tr.IntXhat(1) + IntMarginError;
           Down =  tr.IntXhat(1) - IntMarginError;
           
                                %% DEBUG MODE
                                
                    tr.HitMiss(end + 1) = 0;
                    count = size(tr.HitMiss,2) - 1;
                    
                    tr.TempPMZ(count) = tr.MZP(1,1);
                    tr.TempPInt(count) = tr.IntP(1,1);
    
                    tr.MZConf(count, :) = [Left Right];
                    tr.IntConf(count, :) = [Up Down];
                    
                    tr.MZPred(count) = tr.MZXhat(1);
                    tr.IntPred(count) = tr.IntXhat(1);
                    
                
                    if (tr.ID(1) == 53 && tr.ID(2) == 3) 
                       display('interest'); 
                    end
                    
           %% 
           
          %See what potential centroids maybe in its current scan      
           NarrowIndex = intersect(find(CurrentScan(:,1)>= Left), find(CurrentScan(:,1) <= Right));
           
           
           %just skip it if the current prediction is not in the scan for
           %m/z column => it also is not in the int column.
           
           %the reason for the multiplication is sometimes it returns
           %0 X 1 or sometimes a 1 X 0. So to avoid a false acceptance do
           %this
           if size(NarrowIndex,1)*size(NarrowIndex,2) == 0
               BestPointIdx = -1; %A return value to say that no point was found
               LeastError = -1; %A dummy return value;
               return;
           end
           
           DataPointIdx = ones(size(CurrentScan(NarrowIndex, 1), 1))*Inf;
           Distance = DataPointIdx;
           %evaluate which in the initial pass work
           for i = 1:size(CurrentScan(NarrowIndex), 1)
               
               %intensity check
               if  (Down  <   CurrentScan(NarrowIndex(i),2) && ...
                    Up    >   CurrentScan(NarrowIndex(i),2))
                
                    
                    %the distance measure for choosing the closest point
                    %in the case that more than one point is claimed.
                    DMZ = sqrt((CurrentScan(NarrowIndex(i),1) - ...
                                tr.MZXhat(1))^2)/sqrt(tr.MZP(1,1));
                    DInt = sqrt((CurrentScan(NarrowIndex(i),2) - ....
                                 tr.IntXhat(1))^2)/sqrt(tr.IntP(1,1));
                    Distance(i, 1) = DMZ + DInt;
                    %store the claimed index value
                    DataPointIdx(i, 1) = NarrowIndex(i);
                   
                end      
           end

           %remove preallocated elements that serve as placeholders
           PotentialClaims = find(DataPointIdx ~= Inf);

           %Choose only the best/closest estimate in the 
           %case that there is not one-to-one map     
           %this should happen within Tracker Class because it is
            if size(PotentialClaims,1) > 0      
                ProcessedDistance = Distance(PotentialClaims, 1);
                LeastError = min(ProcessedDistance);
                BestPointIdx = find(Distance == LeastError);
                %choose the first index as a workaround if there is not a 
                %unique minimum least error
                BestPointIdx = NarrowIndex(BestPointIdx(1));
            else
                BestPointIdx = -1; %A return value to say that no point was found
                LeastError = -1; %A dummy return value;
            end   
      end
      function innovateCentroid(tr,ScanNum, DataPointNum, Centroid)
          %@ params: tr->the respective Tracker trect
          %          Centroid->[MZmeasurement IntMeasurement]
          %          ScanNum-> The Current Scan you are in
          %          DataPointNum->the sequential order in which
          %                        Data lies from the current Scan
          %
          % returns: Nothing, function only alters internal model.
          
          %innovation mass
          tr.MZP = inv(inv(tr.MZP) + tr.MZH'*(tr.MZR\tr.MZH));
          tr.MZXhat = tr.MZXhat - tr.MZP*(tr.MZH'/tr.MZR)...
                     *(tr.MZH*tr.MZXhat - Centroid(1));  
           
          %innovation intensity
          tr.IntP = inv(inv(tr.IntP) + tr.IntH'*(tr.IntR\tr.IntH));
          tr.IntXhat = tr.IntXhat - tr.IntP*(tr.IntH'/tr.IntR)...
                     *(tr.IntH*tr.IntXhat - Centroid(2));    
                 
          %STORE DATAPOINT TRACKED
          tr.ScanNumbers(tr.TrackerLength) = ScanNum;
          tr.DataPointNumbers(tr.TrackerLength) = DataPointNum;
        
          %% DEBUG MODE for confidence evolution figure
          tr.MZInno(tr.TrackerLength) = tr.MZXhat(1);
          tr.IntInno(tr.TrackerLength) = tr.IntXhat(1);
          count = size(tr.HitMiss,2) - 1;
          tr.InnoPMZ(count) = tr.MZP(1,1);
          tr.InnoPInt(count) = tr.IntP(1,1);
          tr.HitMiss(end) = 1;
      end
   end %methods
end % classdef

=end






















