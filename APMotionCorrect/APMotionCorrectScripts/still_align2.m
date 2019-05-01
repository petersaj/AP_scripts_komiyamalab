function [stilloffsets]=still_align2(stillimages,maxoffset)
%takes a sequence of STILLIMAGES (numstils x N x N) and returns a matrix
%containing the relative offsets between these images in X and Y based upon
%the shift which gives the maximal cross correlational considered only over
%a range (MAXOFFSET) of offsets around zero offset
  
%pull out the number of still and the size of the image
  numstills=size(stillimages,1);
  N=size(stillimages,2);
  M=size(stillimages,3);
  
  %how many offsets are you going to consider
  Nx=(4*maxoffset)+1;
  Ny=(4*maxoffset)+1;
  
  %we are aligning them all to the first image, so the offset there is
  %zero,zero by definition
  stilloffsets(1,:)=[0 0];
  %it is the reference image
  refimage=squeeze(stillimages(1,:,:));


  %loop over all the other offsets
  for i=2:numstills
 % for i=2:100  
      %pick out the image we are aligning
    testimage=squeeze(stillimages(i,:,:));
    testimage2=interp2(1:2:2*size(refimage,2),1:2:2*size(refimage,1),double(testimage),1:2*size(refimage,2),(1:2*size(refimage,1))','spline',0);
    %initialize the correlation vector
    correlation=zeros(Nx,Ny);
    
    %loop over all shifts
    for shiftx=-2*maxoffset:2*maxoffset
      for shifty=-2*maxoffset:2*maxoffset
          
        %create a hashed index for this combination of x and y offset
        shifthash=((shiftx+2*maxoffset)*Nx) + (shifty+2*maxoffset+1);
        %disp(shifthash);
        
        %cut out the section of reference image for this offset
        %refimage_cut=refimage(1+2*maxoffset:N-maxoffset,1+maxoffset:M- ...
        %                      maxoffset);
           
        testimage_cut=testimage2(1+2*maxoffset+shifty:2:2*N-2*maxoffset+shifty-2,1+2*maxoffset+shiftx:2:2*M-2*maxoffset+shiftx-2);
        testimage_cut=round(testimage_cut);                 
        %cut out the section of test image for this offset
        refimage_cut=refimage(1+maxoffset:N-maxoffset-1,1+maxoffset:M-maxoffset-1);
        
       % testimage_cut=testimage(1+maxoffset+shifty:N-maxoffset+shifty,1+ ...
        %                       maxoffset+shiftx:M-maxoffset+shiftx);
                           
        %subtract out the mean from both
        refimage_cut=refimage_cut-mean(mean(refimage_cut));
        testimage_cut=testimage_cut-mean(mean(testimage_cut));
        
       
        %calculate the cross correlation and store it in the vector using
        %the hash as an index
        correlation(1+2*maxoffset+shifty,1+2*maxoffset+shiftx)=sum(sum(refimage_cut.*testimage_cut));
  
      
        %I considered using least square difference
        %correlation(shifthash)=mean(mean(abs(refimage_cut-testimage_cut)));
      end
    end
         
%          figure(10);
%          imagesc(correlation);
%          pause(.05);
    
    %find the maximum correlation
    [junks,maxx]=max(correlation);
    [junk,maxy]=max(junks);
    [shiftsx,shiftsy]=meshgrid(-2*maxoffset:2*maxoffset);
    
    %reverse the hash to get out the offsets
%     minyshift=mod(minshift_hash,Nx);
%     minxshift=((minshift_hash-minyshift)/Nx);
%     minyshift=minyshift-2*maxoffset-1;
%     minxshift=minxshift-2*maxoffset;
    %save the offsets
    stilloffsets(i,:)=[shiftsy(maxy,maxx(maxy)) shiftsx(maxy,maxx(maxy))];
 
  end
  
  