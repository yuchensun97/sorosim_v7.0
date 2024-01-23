classdef TestFwdKinematics < matlab.unittest.TestCase

    properties
        VLinks = SorosimLink.empty(8, 0); % input links
        Vq = cell(8); % input q
        Vg_true = cell(8); % true g output
        Vrho_true = cell(8); % true rho output
    end
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
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
                                       0 1 1 2 1 0]';
            testCase.VLinks(8).B_rho = [1 2];

            % load test data and ground true results
            d = load("tests\FwdKinematicsData.mat");
            for i=1:4
                testCase.Vq{i} = d.q_in{i};
                testCase.Vg_true{i} = d.g_out{i};
                testCase.Vrho_true{i} = ones(12, 1);
            end
            for i=5:7
                testCase.Vq{i} = d.q_in{1};
                testCase.Vg_true{i} = d.g_out{1};
            end
            testCase.Vq{8} = d.q_in{4};
            testCase.Vg_true{8} = d.g_out{4};
        end
    end
    
    methods(Test)
        % Test methods
        function test_classic_constant(testCase)
            L = SorosimLinkage(testCase.VLinks(1));
            q = testCase.Vq{1};
            g_true = testCase.Vg_true{1};
            rho_true = testCase.Vrho_true{1};
            q_rho = rand;
            [g, rho] = L.FwdKinematics(q, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(rho, rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_linear(testCase)
            L = SorosimLinkage(testCase.VLinks(2));
            q = testCase.Vq{2};
            g_true = testCase.Vg_true{2};
            rho_true = testCase.Vrho_true{2};
            q_rho = rand;
            [g, rho] = L.FwdKinematics(q, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(rho, rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_quardratic(testCase)
            L = SorosimLinkage(testCase.VLinks(3));
            q = testCase.Vq{3};
            q_rho = rand;
            g_true = testCase.Vg_true{3};
            rho_true = testCase.Vrho_true{3};
            [g, rho] = L.FwdKinematics(q, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(rho, rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_combine(testCase)
            L = SorosimLinkage(testCase.VLinks(4));
            q = testCase.Vq{4};
            q_rho = rand;
            g_true = testCase.Vg_true{4};
            rho_true = testCase.Vrho_true{4};
            [g, rho] = L.FwdKinematics(q, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(rho, rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_constant(testCase)
            L = SorosimLinkage(testCase.VLinks(5));
            q_xi = testCase.Vq{5};
            q_rho = -rand;
            g_true = testCase.Vg_true{5};
            rho_true = ones(12, 1) + q_rho;
            rho_true(1) = 1;
            [g, rho] = L.FwdKinematics(q_xi, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(rho, rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_linear(testCase)
            L = SorosimLinkage(testCase.VLinks(6));
            q_xi = testCase.Vq{6};
            q_rho = -rand(2, 1);
            g_true = testCase.Vg_true{6};
            [g, ~] = L.FwdKinematics(q_xi, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
        end

        function test_rho_quardratic(testCase)
            L = SorosimLinkage(testCase.VLinks(7));
            q_xi = testCase.Vq{7};
            q_rho = -rand(3, 1);
            g_true = testCase.Vg_true{7};
            [g, ~] = L.FwdKinematics(q_xi, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
        end

        function test_rho_combine(testCase)
            L = SorosimLinkage(testCase.VLinks(8));
            q_xi = testCase.Vq{8};
            q_rho = -rand(3, 1);
            g_true = testCase.Vg_true{8};
            [g, ~] = L.FwdKinematics(q_xi, q_rho);
            testCase.verifyEqual(g, g_true, 'AbsTol', 1e-10);
        end
    end
    
end
