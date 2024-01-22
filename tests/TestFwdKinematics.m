classdef TestFwdKinematics < matlab.unittest.TestCase

    properties
        VLinks = SorosimLink.empty(8, 0); % input links
        Vg_true = cell(8); % true g output
        Vrho_true = cell(8); % true rho output
    end
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
            rng(0);
            % Create test links
            baseLink = SorosimLink();
            baseLink.r_base = 0.02;
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

            testCase.VLinks(5) = testCase.VLinks(1);
            testCase.VLinks(5).B_rho = [1 0];

            testCase.VLinks(6) = testCase.VLinks(1);
            testCase.VLinks(6).B_rho = [1 1];

            testCase.VLinks(7) = testCase.VLinks(1);
            testCase.VLinks(7).B_rho = [1 2];

            testCase.VLinks(8) = testCase.VLinks(1);
            testCase.VLinks(8).B_xi = [1 1 1 1 1 1;
                                       1 1 1 1 1 1]';
            testCase.VLinks(8).B_rho = [1 2];
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

        function test_rho_constant(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_rho_linear(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_rho_quardratic(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function test_rho_combine(testCase)
            testCase.verifyFail("Unimplemented test");
        end
    end
    
end
