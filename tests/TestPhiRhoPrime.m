classdef TestPhiRhoPrime < matlab.unittest.TestCase

    methods(TestClassSetup)
        % Test class setup
        function setup(testCase)
            rng(34);
        end
    end
    
    methods(Test)
        % Test methods
        function testConstant(testCase)
            n = 0;
            x_in = rand(); % [0, 1]
            x = 2 * x_in - 1;
            prime_true = testCase.legendrePrime(n, x);
            prime = Phi_Prime_Rho_LegendrePolynomial(x_in, 1, n);
            testCase.verifyEqual(prime, prime_true, 'AbsTol', 1e-10); 
        end

        function testLinear(testCase)
            n = 1;
            x_in = rand();
            x = 2 * x_in - 1;
            prime_true = testCase.legendrePrime(n, x);
            prime = Phi_Prime_Rho_LegendrePolynomial(x_in, 1, n);
            testCase.verifyEqual(prime, prime_true, 'AbsTol', 1e-10);
        end

        function testQuadratic(testCase)
            n = 2;
            x_in = rand();
            x = 2 * x_in - 1;
            prime_true = testCase.legendrePrime(n, x);
            prime = Phi_Prime_Rho_LegendrePolynomial(x_in, 1, n);
            testCase.verifyEqual(prime, prime_true, 'AbsTol', 1e-10);
        end

        function testCubic(testCase)
            n = 3;
            x_in = rand();
            x = 2 * x_in - 1;
            prime_true = testCase.legendrePrime(n, x);
            prime = Phi_Prime_Rho_LegendrePolynomial(x_in, 1, n);
            testCase.verifyEqual(prime, prime_true, 'AbsTol', 1e-10);
        end

        function testRandomOrder(testCase)
            n = randi(10);
            x_in = rand();
            x = 2 * x_in - 1;
            prime_true = testCase.legendrePrime(n, x);
            prime = Phi_Prime_Rho_LegendrePolynomial(x_in, 1, n);
            testCase.verifyEqual(prime, prime_true, 'AbsTol', 1e-10);
        end
    end

    methods
        function pn_prime = legendrePrime(testCase,n, x)
        % Legendre polynomial derivative
        % P'_n(x) = nP_{n-1}(x) + xP'_{n-1}(x)
        % returns: [P'_0(x), P'_1(x), ... P'_n(x)]
            pn_prime = zeros(1, n+1);
            pn = zeros(1, n+1);
            if n == 0
                return
            end
            if n == 1
                pn_prime(1) = 0;
                pn_prime(2) = 1;
                return
            end
            pn(1) = 1;
            pn(2) = x;
            pn_prime(1) = 0;
            pn_prime(2) = 1;
            for i = 2:n
                pn(i+1) = ((2*i-1) * x * pn(i) - (i-1) * pn(i-1)) / i;
                pn_prime(i+1) = i * pn(i) + x * pn_prime(i);
            end
        end
    end

end
