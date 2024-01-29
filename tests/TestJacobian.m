classdef TestJacobian < matlab.unittest.TestCase

    properties
        VLinks = SorosimLink.empty(8, 0);
        Vq = cell(1,8);
        VJ_xi_true = cell(1,8);
        Vrho_true = cell(1,8);
    end
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
            % create test links
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

            testCase.VLinks(5) = testCase.VLinks(1);
            testCase.VLinks(5).B_rho = [1 0];

            testCase.VLinks(6) = testCase.VLinks(1);
            testCase.VLinks(6).B_rho = [1 1];

            testCase.VLinks(7) = testCase.VLinks(1);
            testCase.VLinks(7).B_rho = [1 2];

            testCase.VLinks(8) = testCase.VLinks(4);
            testCase.VLinks(8).B_rho = [1 2];

            d = load("tests\JacobianData.mat");
            for ii=1:4
                testCase.Vq{ii} = d.qxi{ii};
                testCase.VJ_xi_true{ii} = d.Jxi{ii};
                testCase.Vrho_true{ii} = zeros(12, 0);
            end
            for ii=5:7
                testCase.Vq{ii} = d.qxi{1};
                testCase.VJ_xi_true{ii} = d.Jxi{1};
            end
            testCase.Vq{8} = d.qxi{4};
            testCase.VJ_xi_true{8} = d.Jxi{4};
        end
    end
    
    methods(Test)
        % Test methods
        function test_classic_constant(testCase)
            L = SorosimLinkage(testCase.VLinks(1));
            q = testCase.Vq{1};
            J_xi_true = testCase.VJ_xi_true{1};
            J_rho_true = testCase.Vrho_true{1};
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_linear(testCase)
            L = SorosimLinkage(testCase.VLinks(2));
            q = testCase.Vq{2};
            J_xi_true = testCase.VJ_xi_true{2};
            J_rho_true = testCase.Vrho_true{2};
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_quardratic(testCase)
            L = SorosimLinkage(testCase.VLinks(3));
            q = testCase.Vq{3};
            J_xi_true = testCase.VJ_xi_true{3};
            J_rho_true = testCase.Vrho_true{3};
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_classic_mixed(testCase)
            L = SorosimLinkage(testCase.VLinks(4));
            q = testCase.Vq{4};
            J_xi_true = testCase.VJ_xi_true{4};
            J_rho_true = testCase.Vrho_true{4};
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_constant(testCase)
            L = SorosimLinkage(testCase.VLinks(5));
            q = testCase.Vq{5};
            J_xi_true = testCase.VJ_xi_true{5};
            J_rho_true = L.Twists(2).B_rho;
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_linear(testCase)
            L = SorosimLinkage(testCase.VLinks(6));
            q = testCase.Vq{6};
            J_xi_true = testCase.VJ_xi_true{6};
            J_rho_true = L.Twists(2).B_rho;
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_quardratic(testCase)
            L = SorosimLinkage(testCase.VLinks(7));
            q = testCase.Vq{7};
            J_xi_true = testCase.VJ_xi_true{7};
            J_rho_true = L.Twists(2).B_rho;
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end

        function test_rho_mixed(testCase)
            L = SorosimLinkage(testCase.VLinks(8));
            q = testCase.Vq{8};
            J_xi_true = testCase.VJ_xi_true{8};
            J_rho_true = L.Twists(2).B_rho;
            [J_xi, J_rho] = L.Jacobian(q);
            testCase.verifyEqual(J_xi, J_xi_true, 'AbsTol', 1e-10);
            testCase.verifyEqual(J_rho, J_rho_true, 'AbsTol', 1e-10);
        end
    end
    
end
