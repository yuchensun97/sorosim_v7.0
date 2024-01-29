classdef TestJacobiandot < matlab.unittest.TestCase

    properties
        VLinks = SorosimLink.empty(4, 0);
        Vq = cell(1,4);
        Vqdot = cell(1,4);
        VJd_true = cell(1,4);
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
            d = load("tests\JacobiandotData.mat");
            for ii=1:4
                testCase.Vq{ii} = d.qxi{ii};
                testCase.Vqdot{ii} = d.qdxi{ii};
                testCase.VJd_true{ii} = d.Jdxi{ii};
            end
        end
    end
    
    methods(Test)
        % Test methods
        
        function test_classic_constant(testCase)
            L = SorosimLinkage(testCase.VLinks(1));
            q = testCase.Vq{1};
            qdot = testCase.Vqdot{1};
            Jd_true = testCase.VJd_true{1};
            Jd = L.Jacobiandot(q, qdot);
            testCase.verifyEqual(Jd, Jd_true, 'AbsTol', 1e-10);
        end

        function test_classic_linear(testCase)
            L = SorosimLinkage(testCase.VLinks(2));
            q = testCase.Vq{2};
            qdot = testCase.Vqdot{2};
            Jd_true = testCase.VJd_true{2};
            Jd = L.Jacobiandot(q, qdot);
            testCase.verifyEqual(Jd, Jd_true, 'AbsTol', 1e-10);
        end

        function test_classic_quardratic(testCase)
            L = SorosimLinkage(testCase.VLinks(3));
            q = testCase.Vq{3};
            qdot = testCase.Vqdot{3};
            Jd_true = testCase.VJd_true{3};
            Jd = L.Jacobiandot(q, qdot);
            testCase.verifyEqual(Jd, Jd_true, 'AbsTol', 1e-10);
        end

        function test_classic_combine(testCase)
            L = SorosimLinkage(testCase.VLinks(4));
            q = testCase.Vq{4};
            qdot = testCase.Vqdot{4};
            Jd_true = testCase.VJd_true{4};
            Jd = L.Jacobiandot(q, qdot);
            testCase.verifyEqual(Jd, Jd_true, 'AbsTol', 1e-10);
        end

    end
    
end
