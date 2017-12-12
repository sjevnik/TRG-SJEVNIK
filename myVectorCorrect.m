function v = myVectorCorrect(nRows,nColumns,nIterations,gridSpacing,interpolationWindow,dataFile)	      
    H = 7.22; %interpolation parameter    
    v = dataFile; %load('D:\Experimental Data\2016-06-27 Test 1\OutputMichael\data0001.piv');
    noc = nColumns;
    nor = nRows;
    gd = gridSpacing;
    itr = nIterations;
    iwin = interpolationWindow;
    R=size(v,1);  
    
    %%From here below is directly copy and pasted from the Jubayer code.
    %%There is only one small part missing and that relates to a variable
    %%named "spvec" which did nothing but store information about the
    %%corrected vectors. I have replaced the single & with && for logical
    %%and as well as replacing | with || for logical or 
    %%Everything else is as I found it unless I
    %%commented my name next to it.
    %%Sorry for the tabbing. Dennis
        % computing vectors at grid points of zero velocity
        for i=R:-1:noc+2
            if (v(i,3)==0 && v(i,4)==0)
                nbr1=0;
                sumu1=0;
                sumv1=0;
                for b=0:noc:(2*noc);
                    for a=0:2;
                        if (((v(i,10)+0.5)/gd)+((b-noc)/noc)>0 && ((v(i,10)+0.5)/gd)-1+((b-noc)/noc)<=nor);
                            if ((i+(a-1)+(b-noc))<=R && (i+(a-1)+(b-noc))>0);
                                if (v((i+(a-1)+(b-noc)),10)==v((i+(b-noc)),10));
                                    if (v((i+(a-1)+(b-noc)),3)~=0 || v((i+(a-1)+(b-noc)),4)~=0);
                                        if (i+(a-1)+(b-noc))~=i
                                            nbr1=nbr1+1;
                                            sumu1=sumu1+v(i+(a-1)+(b-noc),3);
                                            sumv1=sumv1+v(i+(a-1)+(b-noc),4);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                v(i,3) = sumu1/nbr1;
                v(i,4) = sumv1/nbr1;
            end
        end
        
        %==========================================================================
        % First step in replacing the spurious vectors is to replace the spurious
        % vectors by the median of 8 neighbouring vectors

        tot=0;	% to count total vectors
        for b1=1:itr
            for i=1:R
                %         for i=4071	% for numerical varification
                nanv=isnan(v(i,3));
                if nanv==0  % i.e. the value is a real number
                    %             if (v(i,3)~=NaN || v(i,4)~=NaN)
                    tot=tot+1;
                    %if (v(i-47,11)==0 && v(i-47,12)==0)	% for top layer
                    nbr2=0;
                    %sumu=0; %Dennis
                    %sumv=0; %Dennis
                    for b=0:noc:(2*noc);
                        for a=0:2;
                            if (((v(i,10)+0.5)/gd)+((b-noc)/noc)>0 && ((v(i,10)+0.5)/gd)-1+((b-noc)/noc)<=nor);
                                if ((i+(a-1)+(b-noc))<=R && (i+(a-1)+(b-noc))>0);
                                    if ((i+(b-noc))<=R && (i+(b-noc))>0)
                                        if (v((i+(a-1)+(b-noc)),10)==v((i+(b-noc)),10));
                                            nanchk=isnan(v((i+(a-1)+(b-noc)),3));
                                            if nanchk==0    % i.e. the value is a real number
                                                %if (v((i+(a-1)+(b-noc)),3)~=NaN && v((i+(a-1)+(b-noc)),4)~=NaN);
                                                if (i+(a-1)+(b-noc))~=i
                                                    nbr2=nbr2+1;
                                                    arru(nbr2) = v(i+(a-1)+(b-noc),3);	% writing the neighbouring vectors in an array
                                                    arrv(nbr2) = v(i+(a-1)+(b-noc),4);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    mdu = median(arru);
                    mdv = median(arrv);
                    %                 pause;
                    arru=0;
                    arrv=0;
                    % replacing the vectors if they are greater then the set tolerence
                    % the tolerence is based on the median of eight (or less) neighbours

                    gv(i,3)=v(i,3);
                    gv(i,4)=v(i,4);

                    % Step 1.
                    % computing the magnitude of resultant raw vector and median vector and
                    % then computing the ratio of both. If the ratio is within certain limits
                    % specified by "limt" then check for the angular difference between the
                    % median vector and raw vector. If the angular difference is within the
                    % range specified then no correction of vectors other wise correction is
                    % done as specified by the following procedure.
                    % Limits for checking the spurious vectors is based on the magnitude of
                    % the median of neighbouring vectors. If median >= 8.8 then limit=4.4/median.
                    % If 1.5 > median > 8.8 then limit=0.5; if 0.5 > median > 1.5 limit=1.0 and
                    % if median <= 0.5 then limit=2.0. All the raw vectors outside these limits
                    % are replaced by the median vector

                    magmd=sqrt(mdu^2+mdv^2);
                    magrw=sqrt(v(i,3)^2+v(i,4)^2);

                    dirmd=atan2(mdv,mdu);
                    dirrw=atan2(v(i,4),v(i,3));

                    deldir=abs(dirmd-dirrw);

                    rmag=magrw/magmd;

                    limt=0.25;	% limit for magnitude of resultant velocity
                    thl=15*pi/180;	% limit for angular difference
                    thh=(360-15)*pi/180;	% limit for angular difference

                    mxval=2.0;	% maximum value of limit for intermedate range velocities
                    mxval1=4.0; % maximun value of limit for smaller velocities

                    mvel1=1.168; % (equal to 6.25 cm/s)velocity limit for the critria limit=factor/abs(median)
                    mvel2=0.187; % (equal to 1.0 cm/s)

                    fctr=2.336;	% a constant value to set the limit (equivalent to 12.5)

                    if rmag>=(1-limt) && rmag<=(1+limt)
                        %i; Dennis
                        if deldir<=thl || deldir>=thh
                            %i; Dennis
                        else	% to execute if magnitude is within limits but angular difference is
                            % outside the limits

                            if mdu==0
                                mdu=0.1;
                            end
                            if mdv==0
                                mdv=0.1;
                            end

                            % large velocities
                            if abs(mdu)>=mvel1
                                lmu=fctr/abs(mdu);
                                % applying the tolerence limits
                                if mdu>0
                                    if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);	% computing the x-location for the corrected vector
                                        v(i,2)=v(i,10)+(gv(i,4)/2);	% computing the y-location for the corrected vector
                                    end
                                end
                                if mdu<0
                                    if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end

                            % intermediate velocities
                            if abs(mdu)>mvel2 && abs(mdu)<mvel1
                                lmu=mxval;
                                % applying the tolerence limits
                                if mdu>0
                                    if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                                if mdu<0
                                    if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end

                            % small velocities
                            if abs(mdu)<=mvel2
                                lmu=mxval1;
                                % applying the tolerence limits
                                if mdu>0
                                    if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                                if mdu<0
                                    if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end
                            
                            % V comp.
                            % large velocities
                            if abs(mdv)>=mvel1
                                lmv=fctr/abs(mdv);
                                % applying the tolerence limits
                                if mdv>0
                                    if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                                if mdv<0
                                    if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end
                            
                            % intermediate velocities
                            if abs(mdv)>mvel2 && abs(mdv)<mvel1
                                lmv=mxval;
                                % applying the tolerence limits
                                if mdv>0
                                    if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                                if mdv<0
                                    if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end

                            % small velocities
                            if abs(mdv)<=mvel2
                                lmv=mxval1;
                                % applying the tolerence limits
                                if mdv>0
                                    if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                                if mdv<0
                                    if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                        gv(i,3)=mdu;
                                        gv(i,4)=mdv;
                                        v(i,1)=v(i,9)+(gv(i,3)/2);
                                        v(i,2)=v(i,10)+(gv(i,4)/2);
                                    end
                                end
                            end                          
                        end  % end of if statement for angular difference

                    else	% to execute correction step if magnitude is outside the limits
                        if mdu==0
                            mdu=0.1;
                        end
                        if mdv==0
                            mdv=0.1;
                        end

                        % large velocities
                        if abs(mdu)>=mvel1
                            lmu=fctr/abs(mdu);
                            % applying the tolerence limits
                            if mdu>0
                                if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdu<0
                                if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end

                        % intermediate velocities
                        if abs(mdu)>mvel2 && abs(mdu)<mvel1
                            lmu=mxval;
                            % applying the tolerence limits
                            if mdu>0
                                if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdu<0
                                if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end

                        % small velocities
                        if abs(mdu)<=mvel2
                            lmu=mxval1;
                            % applying the tolerence limits
                            if mdu>0
                                if v(i,3)>(mdu*(1+lmu)) || v(i,3)<(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdu<0
                                if v(i,3)<(mdu*(1+lmu)) || v(i,3)>(mdu*(1-lmu))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end
                        
                        % V comp.
                        % large velocities
                        if abs(mdv)>=mvel1
                            lmv=fctr/abs(mdv);
                            % applying the tolerence limits
                            if mdv>0
                                if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdv<0
                                if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end
                        
                        % intermediate velocities
                        if abs(mdv)>mvel2 && abs(mdv)<mvel1
                            lmv=mxval;
                            % applying the tolerence limits
                            if mdv>0
                                if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdv<0
                                if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end
                        
                        % small velocities
                        if abs(mdv)<=mvel2
                            lmv=mxval1;
                            % applying the tolerence limits
                            if mdv>0
                                if v(i,4)>(mdv*(1+lmv)) || v(i,4)<(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                            if mdv<0
                                if v(i,4)<(mdv*(1+lmv)) || v(i,4)>(mdv*(1-lmv))
                                    gv(i,3)=mdu;
                                    gv(i,4)=mdv;
                                    v(i,1)=v(i,9)+(gv(i,3)/2);
                                    v(i,2)=v(i,10)+(gv(i,4)/2);
                                end
                            end
                        end                       
                    end	 % end of if statement for magnitude ratio
                end	% end for non-zero u or v loop
            end	% end for i loop


            %if b1==itr
            %quiver2(v(:,9),v(:,10),gv(:,3),gv(:,4),s,'b')
            %end
            %if b1==3
            %quiver(v(:,9),v(:,10),gv(:,3),gv(:,4),s,'c')
            %end

            %    cct=0;
            %    for k=R:-1:i
            %        cct=cct+1;
            %        nanv1=isnan(v(cct,3));
            %        if nanv1==0
            %            v(cct,3)=gv(cct,3);
            %            v(cct,4)=gv(cct,4);
            %        end
            %    end

            for k=1:R
                nanv1=isnan(v(k,3));
                if nanv1==0
                    v(k,3)=gv(k,3);
                    v(k,4)=gv(k,4);
                elseif nanv1==1
                    v(k,11:12)=NaN;
                end
            end
        end	% end for iteration
        % pause;
        gv=0;

        %=====================================================================================
        % interpolating the velocity vectors onto the grid points using AGW
        % for i=41;
        for i=1:R
            nanv2=isnan(v(i,3));
            if nanv2==0  % i.e. the value is a real number
                %    if (v(i,3)~=0 || v(i,4)~=0)
                sumai=0;
                sumau=0;
                sumav=0;
                for b=0:noc:(iwin*noc);
                    for a=0:iwin;
                        if (((v(i,10)+0.5)/gd)+((b-((iwin/2)*noc))/noc)>0 && ((v(i,10)+0.5)/gd)-1+((b-((iwin/2)*noc))/noc)<=nor); % data lies with in first and last rows
                            if ((i+(a-(iwin/2))+(b-((iwin/2)*noc)))<=R && (i+(a-(iwin/2))+(b-((iwin/2)*noc)))>0); % lies within 1 to R
                                if ((i+(b-((iwin/2)*noc)))>0 && (i+(b-((iwin/2)*noc)))<=R)
                                    if (v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),10)==v((i+(b-((iwin/2)*noc))),10));	% in same row
                                        nanchk2 = isnan(v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),3));
                                        if nanchk2==0
                                            %                      	if (v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),3)~=0 || v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),4)~=0); % considering only non-zero vectors for averaging.
                                            ai=exp(-abs((v(i,9)-v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),1))^2+(v(i,10)-v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),2))^2)/(H^2)); %(x-xi)^2+(y-yi)^2
                                            sumai=sumai+ai;
                                            au=ai*v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),3);	% ai*Ui
                                            av=ai*v((i+(a-(iwin/2))+(b-((iwin/2)*noc))),4);	% ai*Vi
                                            sumau=sumau+au;
                                            sumav=sumav+av;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                v(i,11)=sumau/sumai;	% U=E(ai*Ui)/E(ai)
                v(i,12)=sumav/sumai;
            end
        end
end