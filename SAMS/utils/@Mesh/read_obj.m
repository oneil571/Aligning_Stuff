  function [V,F,Fs] = read_obj( filename )
% function [V,F] = read_obj( filename )
%
% Loads a mesh in Wavefront OBJ format.

   fp = fopen( filename, 'r' );
   if( fp == -1 )
      disp( sprintf( 'Error: could not read mesh file "%s"\n', filename ));
      return;
   end

   V = [];
   F = [];
   line = fgets(fp); line = fgets(fp); line = fgets(fp); line = fgets(fp);     %skip down to relevant part
   while( ~feof( fp ))
       line = fgets(fp);
       line = strsplit(line);           %edit to fix problems with reading texture coordinates
       
        if strcmp(line{1},'v') % vertex
              V = [V [str2double(line{2}); str2double(line{3}); str2double(line{4})]];
        elseif strcmp(line{1},'f')
           if( line{1}== 'f' ) % face
              for q = 2:4
                  face1 = strsplit(line{2},'/');
                  face1 = str2double(face1{1});
                  face2 = strsplit(line{3},'/');
                  face2 = str2double(face2{1});
                  face3 = strsplit(line{4},'/');
                  face3 = str2double(face3{1});
              end
              F = [F [face1;face2;face3]];
           end
       end
   end
   F = F-min(min(F))+1;
   fclose( fp );
   Fs = [];
end

