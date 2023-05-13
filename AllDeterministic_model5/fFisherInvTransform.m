function [dInvFTr] = fFisherInvTransform(vPar)
	dInvFTr = ((exp(2 * vPar) - 1) ./ (exp(2 * vPar) + 1));
end

