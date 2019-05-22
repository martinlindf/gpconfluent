function C = mtimes3(A, B)

s=max(ndims(A), ndims(B)); 
idxA = 1:s+1; 
idxA([2 end]) = idxA([end 2]); 
idxB = 1:s+1; 
idxB([1 end]) = idxB([end 1]); 
C = sum(bsxfun(@times, permute(A, idxA), permute(B, idxB)), s+1);
