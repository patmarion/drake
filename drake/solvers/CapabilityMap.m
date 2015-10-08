classdef CapabilityMap
  
  properties
    map
    reachabilityIndex
    sphCenters
    sphDiameter
    nSamples
    angTolerance
    posTolerance
    urdf
    nSph
    nPointsPerSphere
    rootLink
    rootPoint
    endEffectorLink
    endEffectorPoint
    baseLink
  end
  
  methods
    
    function obj = CapabilityMap(matFile)
      if nargin > 0
        obj = obj.generateFromFile(matFile);
      end
    end
    
    function obj = generateFromFile(obj, matFile)
      vars = load(matFile);
      obj.map = vars.map;
      obj.reachabilityIndex = vars.reachabilityIndex;
      obj.sphCenters = vars.sphCenters;
      obj.sphDiameter = vars.options.sphDiameter;
      obj.nSamples = vars.options.nSamples;
      obj.angTolerance = vars.options.angTolerance;
      obj.posTolerance = vars.options.posTolerance;
      obj.urdf = vars.options.urdf;
      obj.nSph = size(obj.map, 1);
      obj.nPointsPerSphere = size(obj.map, 2);
      obj.rootLink = vars.options.rootLink;
      obj.rootPoint = vars.options.rootPoint;
      obj.endEffectorLink = vars.options.endEffectorLink;
      obj.endEffectorPoint = vars.options.endEffectorPoint;
      obj.baseLink = vars.options.baseLink;
    end
    
    function points = findPointsFromDirection(obj, direction, threshold)
      [P, frames] = distributePointsOnSphere(obj.nPointsPerSphere);
      points = false(obj.nPointsPerSphere, 1);
      sphere();
      hold on
      plot3(P(1,:), P(2,:), P(3,:), 'r.')
      for p = 1:obj.nPointsPerSphere
        if acos(frames(p,:,3)*direction)/norm(direction) <= threshold
          points(p) = true;
%           plot3([P(1,p), P(1,p) - frames(p, 1, 3)'], [P(2,p), P(2,p) - frames(p, 2, 3)'], [P(3,p), P(3,p) - frames(p, 3, 3)'], 'b')
        end
      end
    end

    function [spheres, obj] = findSpheresFromDirection(obj, direction, minSph, maxSph, sagittalAngle,...
        transverseAngle, sagittalWeight, transverseWeight)
      if nargin > 4
        obj = obj.prune(sagittalAngle, transverseAngle, sagittalWeight, transverseWeight, 0);
      end
      spheres = false(obj.nSph, 1);
      ns = 0;
      threshold = 0;
      while ns < minSph
        points = obj.findPointsFromDirection(direction, threshold);
        for s = 1:obj.nSph
          if any(obj.map(s, points))
            spheres(s) = true;
          end
        end
        ns = nnz(spheres);
        threshold = threshold + pi/50;
      end
      if ns > maxSph
        reachabilityWeight = 0.5;
        while ns > maxSph
          obj = obj.prune(sagittalAngle, transverseAngle, sagittalWeight, transverseWeight, reachabilityWeight);
          reachabilityWeight = reachabilityWeight + 0.5;
          points = obj.findPointsFromDirection(direction, threshold - pi/50);
          spheres = false(obj.nSph, 1);
          for s = 1:obj.nSph
            if any(obj.map(s, points))
              spheres(s) = true;
            end
          end
          ns = nnz(spheres);
        end
      end
    end
    
    function drawCapabilityMap(obj, direction, minSph, maxSph)
      lcmClient = LCMGLClient('CapabilityMap');
      [spheres, obj] = obj.findSpheresFromDirection(direction, minSph, maxSph, 0, 0, 2, 1.5);
      for sph = 1:obj.nSph
        if spheres(sph)
          lcmClient.sphere(obj.sphCenters(:,sph), obj.sphDiameter/2, 20, 20);
        end
      end
      disp(nnz(spheres))
      lcmClient.switchBuffers();
    end
    
    function obj = prune(obj, sagittalAngle,...
        transverseAngle, sagittalWeight, transverseWeight, reachabilityWeight)
      
      Dmax = max(obj.reachabilityIndex);
      indices = [];
      
      for sph = 1:obj.nSph
        sa = atan2(obj.sphCenters(3,sph), obj.sphCenters(1,sph));
        ta = atan2(obj.sphCenters(2,sph), obj.sphCenters(1,sph));
        sagittalCost = sagittalWeight * abs(sa - sagittalAngle);
        transverseCost = transverseWeight * abs(ta - transverseAngle);
        reachabilityCost = reachabilityWeight * (Dmax - obj.reachabilityIndex(sph));
        if sqrt(sagittalCost^2 + transverseCost^2) + reachabilityCost < 2
          indices(end + 1) = sph;
        end
      end
      obj.nSph = length(indices);
      obj.reachabilityIndex = obj.reachabilityIndex(indices);
      obj.map = obj.map(indices, :);
      obj.sphCenters = obj.sphCenters(:, indices);
    end
    
  end
    
  methods (Static)
    
    function [P, frames] = distributePointsOnSphere(N)
      k = 1:N;
      h = -1 + 2*(k-1)/(N-1);
      theta = acos(h);
      phi = zeros(1,N);
      for i = 2:N-1
        phi(i) = mod(phi(i-1) + 3.6/(sqrt(N) * sqrt(1 - h(i)^2)), 2*pi);
      end
      x = cos(phi).*sin(theta);
      y = sin(phi).*sin(theta);
      z = cos(theta);
      P = [x; y; z];
      frames = zeros(N, 3, 3);
      for p = 1:N
        frame = zeros(3);
        frame(1:3,3) = -P(:,p);
        if abs(frame(1,3)) <= 1e-10 && abs(frame(2,3)) <= 1e-10
          frame(:,1) = [sign(frame(3,3)); 0; 0];
        else
          frame(2,1) = sqrt(frame(1,3)^2/(frame(2,3)^2 + frame(1,3)^2));
          frame(1,1) = -sign(frame(1,3)*frame(2,3))*sqrt(1-frame(2,1)^2);
        end
        frame(:,2) = cross(frame(:,3), frame(:,1));
        frames(p, :, :) = frame;
      end
    end
    
  end
  
end