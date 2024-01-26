classdef TestPhiRhoPrime < matlab.unittest.TestCase

    properties
        x = linspace(-1, 1, 100);
    end
    
    methods(TestMethodSetup)
        % Setup for each test
        function pn_prime = legrendrePrime(n, x)
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
            for i = 2:n+1
                pn(i+1) = ((2*i-1) * x * pn(i) - (i-1) * pn(i-1)) / i;
                pn_prime(i+1) = i * pn(i) + x * pn_prime(i);
            end
        end

        function pn_prime = legrendrePrimeFunc(n)
           pn_prime = @(x) legrendrePrime(n, x);
        end

    end
    
    methods(Test)
        % Test methods
        function testConstant(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function testLinear(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function testQuadratic(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function testCubic(testCase)
            testCase.verifyFail("Unimplemented test");
        end

        function testRandomOrder(testCase)
            testCase.verifyFail("Unimplemented test");
        end
    end
end
