function out_im = cropCenter(in_im, crop_size)
% Input: Image, size of crop 
% Output: cropped image
[m,n,c] = size(in_im);
m_h = m/2;
n_h = n/2;
c_h = crop_size/2;

out_im = in_im(m_h-c_h(1)+1:m_h+c_h(1),n_h-c_h(2)+1:n_h+c_h(2),:);
end