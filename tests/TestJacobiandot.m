classdef TestJacobiandot < matlab.unittest.TestCase

    properties
        VLinks = SorosimLink.empty(8, 0);
        Vq = cell(8);
        VJ_xi_true = cell(8);
    end
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
            %create test links
            baseLink = SorosimLink();
            testCase.VLinks(1) = baseLink;

            testCase.VLinks(2) = testCase.VLinks(1);
            testCase.VLinks(2).B_xi = [1 1 1 1 1 1;
                                       1 1 1 1 1 1]';
            
            testCase.VLinks(3) = testCase.VLinks(1);
            testCase.VLinks(3).B_xi = [1 1 1 1 1 1;
                                       2 2 2 2 2 2]';
            
            testCase.VLinks(4) = testCase.VLinks(1);
            testCase.VLinks(4).B_xi = [1 1 1 1 1 1;
                                       0 1 1 2 1 0]';
        end
    end
    
    methods(Test)
        % Test methods
        
        function test_classic_constant(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_classic_linear(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_classic_quardratic(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_classic_combine(testCase)
            testCase.verifyFail("Unimplemented test");
        end

    end
    
end
