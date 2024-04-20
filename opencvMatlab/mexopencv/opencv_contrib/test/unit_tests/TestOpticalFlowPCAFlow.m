classdef TestOpticalFlowPCAFlow
    %TestOpticalFlowPCAFlow

    methods (Static)
        function test_1
            %TODO: crashes MATLAB sometimes!
            if true
                error('mexopencv:testskip', 'todo');
            end

            im1 = 255*uint8([...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 1 1 1 0 0 0 0;...
                0 0 0 1 0 1 0 0 0 0;...
                0 0 0 1 1 1 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
                0 0 0 0 0 0 0 0 0 0;...
            ]);
            im2 = circshift(im1, [0 1]);
            alg = cv.OpticalFlowPCAFlow();
            flow = alg.calc(im1, im2);
            validateattributes(flow, {'single'}, ...
                {'3d', 'size',[size(im1,1) size(im1,2) 2]});
        end

        function test_2
            prevImg = cv.imread(fullfile(mexopencv.root(),'test','RubberWhale1.png'), ...
                'Grayscale',true, 'ReduceScale',2);
            nextImg = cv.imread(fullfile(mexopencv.root(),'test','RubberWhale2.png'), ...
                'Grayscale',true, 'ReduceScale',2);
            alg = cv.OpticalFlowPCAFlow();
            flow = alg.calc(prevImg, nextImg);
            validateattributes(flow, {'single'}, ...
                {'3d', 'size',[size(prevImg,1) size(prevImg,2) 2]});
        end
    end

end
