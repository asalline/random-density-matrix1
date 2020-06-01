%Random density matrix generator for qubits

%This script generates random density matrices for pure states of qubits.

%Author: Antti Sällinen
%Last update: 28.5.2020

%First script generates random state vectors correlating to amount of
%wanted qubits.

clear, clc

qubits = input("How many qubits do you want to use?: ")
rank = input("What rank of density matrix you want to use?: ")
Y=1;, y=Y;
N=0;, n=N;
real = input("Do you want complex valued density matrix? (Y/N): ");

%Next this script generates random state vectors for each qubit. The base
%is (1 0)' and (0 1)'. For probability amplitudes there is the requirement
%that |alpha|^2 + |beta|^2 = 1.

ket_0 = [1,0]';
ket_1 = [0,1]';

%"state" contains all the generated qubits. One can look at particular
%qubit by command "state{i}" where i = 1, 2,...,n.

state = cell(1, qubits);
ranks = cell(1, rank);
mixed = cell(1, qubits);

%Script generates either real or complex valued state vectors. Complex
%valued ones have randomly generated signs of the imaginary part.

if rank == 1
    if real == 0
        for k = 1:qubits
            alpha = rand;
            beta = sqrt(1-alpha^2);
            state{k} = alpha * ket_0 + beta * ket_1;
        end
    else
        for k = 1:qubits
            sign = round(rand);
            if sign == 1
                alpha = rand - 1i*rand;
                beta = sqrt(1-alpha^2);
                state{k} = alpha * ket_0 + beta * ket_1;
            else
                alpha = rand + 1i*rand;
                beta = sqrt(1-alpha^2);
                state{k} = alpha * ket_0 + beta * ket_1;
            end
        end
    end
else
    if real == 0
        for k = 1:qubits
            alpha = rand;
            beta = sqrt(1-alpha^2);
            state{k} = alpha * ket_0 + beta * ket_1;
        end
        for l = 1:rank-1
            alpha = rand;
            beta = sqrt(1-alpha^2);
            ranks{l} = alpha * ket_0 + beta * ket_1;
            mixed{l} = kron(state{l},ranks{l})
        end
    else 
        for k = 1:qubits
            sign = round(rand);
            if sign == 1
                alpha = rand - 1i*rand;
                beta = sqrt(1-alpha^2);
                state{k} = alpha * ket_0 + beta * ket_1;
            else
                alpha = rand + 1i*rand;
                beta = sqrt(1-alpha^2);
                state{k} = alpha * ket_0 + beta * ket_1;
            end
        end
    end

%For pure states this part of the script uses tensor product to form the
%density matrix of generated random states. The dimensions of this density
%matrix are 2^n by 2^n, where n is the amount of qubits used.
%After the for-loop the density matrix is normalised to unit by its trace.

if qubits == 1
    rho = state{1} * state{1}'
else
    if rank ~= 1
        prob = rand(1, qubits);
        prob = prob / sum(prob);
        rho = zeros(2^qubits, 2^qubits);
        for k = 1:qubits
            mixed{k}
        end
    else
        for j = 1:(length(state)-1)
            if j+1 == 2
            pure = state{j} * state{j+1}';
            pure = reshape(pure', [] , 1);
            else
                pure = pure * state{j+1}';
                pure = reshape(pure', [], 1);
            end
        end
    end
end
    rho = pure * pure';
    Trace = trace(rho);
    rho = rho/Trace
    Trace = trace(rho)
end

%Below this is some functions that helped to check if this script works
%properly. Because of density matrix the trace should be equals to unit and
%all the eigenvalues should be non negative.

%Trace = trace(rho);
%eigenvalues = eig(rho);
