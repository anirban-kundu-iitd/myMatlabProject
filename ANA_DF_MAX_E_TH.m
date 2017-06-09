clc;
clear all;

snr_dB=0:2:40;
no_simu=10^6;
gamma_th_dB=3;
gamma_th=10^(gamma_th_dB/10);

main_SNR0_dB=3;
main_SNR0=10^(main_SNR0_dB/10);
beta0=1/main_SNR0;

evsd_SNR0_dB=-3;
evsd_SNR0=10^(evsd_SNR0_dB/10);
alpha0=1/evsd_SNR0;

%%-------------------- alpha same-----------------------------
% evsd_SNR1_dB=[3];
% evsd_SNR1_dB=[3 3];
% evsd_SNR1_dB=[3 3 3];
% evsd_SNR1_dB=[3 3 3 3];
%-------------------- alpha different------------------------
evsd_SNR1_dB=[6 9];
% evsd_SNR1_dB=[0 3 6 9];
% %     -----------------------------------------------------

evsd_SNR1=10.^(evsd_SNR1_dB/10);
no_relay=length(evsd_SNR1_dB);
alpha1(1:no_relay)=1./evsd_SNR1;


R=1;
R_s=2^(2*R);
i_snr=1;

% % % %  ----------------beta1 fixed ----------------
% beta1_SNR_dB=30;
% beta1_SNR=10^(beta1_SNR_dB/10);
% beta1(1:no_relay)=1/beta1_SNR;
% %     ---------------------------------------------------

relay_set=1:no_relay;
for total_SNR_dB=snr_dB

    % ----------------beta balanced-----------------------------
    %     beta1_ratio=0.5 ;
    %     beta2_ratio=beta1_ratio;
    %
    %     beta1(1:no_relay)=1./(beta1_ratio.*10^(total_SNR_dB/10));
    %     beta2(1:no_relay)=1./(beta2_ratio.*10^(total_SNR_dB/10));

    % ----------------beta different-----------------------------
    beta1_ratio=[0.20 .30];
%         beta1_ratio=[0.05 0.1 0.15 0.20];
    beta2_ratio=beta1_ratio;

    beta1(1:no_relay)=1./(beta1_ratio.*10^(total_SNR_dB/10));
    beta2(1:no_relay)=1./(beta2_ratio.*10^(total_SNR_dB/10));
    % %     ------------------------------------------------------



    % ---------------- beta1 fixed ----------------
    %         beta2_ratio=[1];
    %         beta2_ratio=[1 1];
    %         beta2_ratio=[1 1 1];
    %         beta2_ratio=[1 1 1 1];
    %
    %         beta2(1:no_relay)=1./(beta2_ratio.*10^(total_SNR_dB/10));
    %     ------------------------------
    %     total_SNR=10.^(total_SNR_dB/10);
    %     beta2=[ beta2 beta0];    %required; diff from simu be careful
    %     final_outg=0;

    beta_D=beta2;
    alpha_E=alpha1;
    % ---------------------------------


    %     rel_sel_temp=1;
    beta=alpha_E;
    temp_out_sel=0;
    temp_asymp_out=0;
    final_outg=0;
         temp_out=0;
    prob_nosel=1;

     %     %     for no_relay=1
%     if no_relay==1
%       prob_sel=exp(-beta1*gamma_th);
       
%         temp_out_sel=1-beta0.*alpha0.*alpha1*exp(-beta2*(R_s-1))...
%             /((beta0-beta2).*(alpha0+R_s*beta2).*(alpha1+R_s*beta2))...
%             -beta2.*alpha0.*alpha1*exp(-beta0*(R_s-1))...
%             /((beta2-beta0).*(alpha0+R_s*beta0).*(alpha1+R_s*beta0));
%         sec_outg(i_snr)=prob_sel*temp_out_sel+(1-prob_sel)*(1-alpha0*exp(-(R_s-1)*beta0)/(alpha0+R_s*beta0));

%         ---------------------------------

%     else
    for i_relay=1:no_relay
        if i_relay==1
            continue
        end
        chosen_set=nchoosek(relay_set,i_relay)
        size_chosen_set=nchoosek(no_relay,i_relay)
        append=(no_relay+1)*ones(size_chosen_set,1)
        chosen_set_append=[chosen_set ]%append
%         temp_out=0;
        
        for i_setrow=1:size_chosen_set
            for i_m=1:size_chosen_set-1
                if i_m==1
                    for i_1=1:size_chosen_set
                        if i_1~=i_setrow
                            temp_beta_1=beta(chosen_set_append(i_1));

                            %----------------------------------------------------------------------
                            %     beta_D=beta1+beta2;
                            %     alpha_E=beta1+alpha1;

                            B_1=beta_D(chosen_set_append(i_setrow))*beta0/(beta_D(chosen_set_append(i_setrow))-beta0);
                            B_2=beta_D(chosen_set_append(i_setrow))*beta0/(beta0-beta_D(chosen_set_append(i_setrow)));
                            num_A=B_1*exp(-beta0*(R_s-1));
                            den_A1=alpha_E(chosen_set_append(i_setrow))+R_s*beta0;
                            den_A2=alpha0+R_s*beta0;

                            num_B=B_2*exp(-beta_D(chosen_set_append(i_setrow))*(R_s-1));
                            den_B1=alpha_E(chosen_set_append(i_setrow ))+R_s*beta_D(chosen_set_append(i_setrow ));
                            den_B2=alpha0+R_s*beta_D(chosen_set_append(i_setrow ));

                            num_C=num_A;
                            den_C1=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1+R_s*beta0;
                            den_C2=alpha0+R_s*beta0;

                            num_D=num_B;
                            den_D1=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1+R_s*beta_D(chosen_set_append(i_setrow ));
                            den_D2=alpha0+R_s*beta_D(chosen_set_append(i_setrow ));

                            A=num_A/(den_A1*den_A2);
                            B=num_B/(den_B1*den_B2);
                            C=num_C/(den_C1*den_C2);
                            D=num_B/(den_D1*den_D2);


                            E=B_1/(beta0*alpha0);
                            F=B_2/(beta_D(i_relay)*alpha0);
                            G=num_A/(beta0*den_A2);
                            H=num_B/(beta_D(i_relay)*den_B2);

                            I_1=R_s*alpha0*(A+B-C-D);
                            I_2=R_s*temp_beta_1*alpha0*(C+D)/(alpha_E(chosen_set_append(i_setrow ))+temp_beta_1);
                            I_3=temp_beta_1*alpha0*(E+F-G-H)/(alpha_E(chosen_set_append(i_setrow ))+temp_beta_1);

                            temp_out_sel=temp_out_sel-(-1)^i_m*(I_1+I_2+I_3);
                            %---------------------------------------------------

                            num_A11=1;
                            den_A11=alpha_E(chosen_set_append(i_setrow ));

                            num_A21=1;
                            den_A21=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;

                            num_A31=alpha0;
                            num_A32=exp(-beta0*(R_s-1));
                            den_A31=alpha0+R_s*beta0;
                            den_A32=alpha_E(chosen_set_append(i_setrow ))+R_s*beta0;

                            num_A41=alpha0;
                            num_A42=exp(-beta0*(R_s-1));
                            den_A41=alpha0+R_s*beta0;
                            den_A42=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1+R_s*beta0;

                            num_B11=alpha_E(chosen_set_append(i_setrow ));
                            den_B11=(alpha_E(chosen_set_append(i_setrow ))+temp_beta_1).^2;

                            num_B21=alpha0;
                            num_B22=alpha_E(chosen_set_append(i_setrow ));
                            num_B23=exp(-beta0*(R_s-1));
                            den_B21=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;
                            den_B22=alpha0+R_s*beta0;
                            den_B23=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1+R_s*beta0;

                            num_C11=R_s;
                            num_C12=temp_beta_1;
                            den_C11=alpha0;
                            den_C12=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;

                            num_C21=temp_beta_1;
                            num_C22=alpha0;
                            num_C23=exp(-beta0*(R_s-1));
                            den_C21=beta0;
                            den_C22=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;
                            den_C23=alpha0+R_s*beta0;

                            num_C31=temp_beta_1;
                            den_C31=beta0;
                            den_C32=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;

                            num_C41=R_s-1;
                            num_C42=temp_beta_1;
                            den_C41=alpha_E(chosen_set_append(i_setrow ))+temp_beta_1;

                            I_1_as=2*R_s*(num_A11/den_A11-num_A21/den_A21...
                                -num_A31*num_A32/(den_A31*den_A32) +num_A41*num_A42/(den_A41*den_A42));
                            I_2_as=2*R_s*(num_B11/den_B11...
                                -num_B21*num_B22*num_B23/(den_B21*den_B22*den_B23));
                            I_3_as=2*(num_C11*num_C12/(den_C11*den_C12)...
                                +num_C21*num_C22*num_C23/(den_C21*den_C22*den_C23)...
                                -num_C31/(den_C31*den_C32)+num_C41*num_C42/den_C41);

                            %                             temp_asymp_out=temp_asymp_out-(-1)^i_m*(I_1_as+I_2_as+I_3_as)/(1/beta1(i_setrow));

                        end  % i_1~=i_setrow
                    end
                end   %if i_m==1


                if i_m==2
                    for i_1=1:size_chosen_set-1
                        for i_2=i_1+1:size_chosen_set
                            if i_1~=i_setrow & i_2~=i_setrow
                                temp_beta_1=beta(chosen_set_append(i_1))+beta(chosen_set_append(i_2));

                                %----------------------------------------------------------------------
                                %     beta_D=beta1+beta2;
                                %     alpha_E=beta1+alpha1;

                                B_1=beta_D(chosen_set_append(i_setrow))*beta0/(beta_D(chosen_set_append(i_setrow))-beta0);
                                B_2=beta_D(chosen_set_append(i_setrow))*beta0/(beta0-beta_D(chosen_set_append(i_setrow)));
                                num_A=B_1*exp(-beta0*(R_s-1));
                                den_A1=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;
                                den_A2=alpha0+R_s*beta0;

                                num_B=B_2*exp(-beta_D(chosen_set_append(i_setrow))*(R_s-1));
                                den_B1=alpha_E(chosen_set_append(i_setrow))+R_s*beta_D(chosen_set_append(i_setrow));
                                den_B2=alpha0+R_s*beta_D(chosen_set_append(i_setrow));

                                num_C=num_A;
                                den_C1=alpha_E(chosen_set_append(i_setrow))+temp_beta_1+R_s*beta0;
                                den_C2=alpha0+R_s*beta0;

                                num_D=num_B;
                                den_D1=alpha_E(chosen_set_append(i_setrow))+temp_beta_1+R_s*beta_D(chosen_set_append(i_setrow));
                                den_D2=alpha0+R_s*beta_D(chosen_set_append(i_setrow));

                                A=num_A/(den_A1*den_A2);
                                B=num_B/(den_B1*den_B2);
                                C=num_C/(den_C1*den_C2);
                                D=num_B/(den_D1*den_D2);


                                E=B_1/(beta0*alpha0);
                                F=B_2/(beta_D(i_relay)*alpha0);
                                G=num_A/(beta0*den_A2);
                                H=num_B/(beta_D(i_relay)*den_B2);

                                I_1=R_s*alpha0*(A+B-C-D);
                                I_2=R_s*temp_beta_1*alpha0*(C+D)/(alpha_E(chosen_set_append(i_setrow))+temp_beta_1);
                                I_3=temp_beta_1*alpha0*(E+F-G-H)/(alpha_E(chosen_set_append(i_setrow))+temp_beta_1);

                                temp_out_sel=temp_out_sel-(-1)^i_m*(I_1+I_2+I_3);
                                %---------------------------------------------------

                                num_A11=1;
                                den_A11=alpha_E(chosen_set_append(i_setrow,i_1));

                                num_A21=1;
                                den_A21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                num_A31=alpha0;
                                num_A32=exp(-beta0*(R_s-1));
                                den_A31=alpha0+R_s*beta0;
                                den_A32=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;

                                num_A41=alpha0;
                                num_A42=exp(-beta0*(R_s-1));
                                den_A41=alpha0+R_s*beta0;
                                den_A42=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                num_B11=alpha_E(chosen_set_append(i_setrow,i_1));
                                den_B11=(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1).^2;

                                num_B21=alpha0;
                                num_B22=alpha_E(chosen_set_append(i_setrow,i_1));
                                num_B23=exp(-beta0*(R_s-1));
                                den_B21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                den_B22=alpha0+R_s*beta0;
                                den_B23=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                num_C11=R_s;
                                num_C12=temp_beta_1;
                                den_C11=alpha0;
                                den_C12=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                num_C21=temp_beta_1;
                                num_C22=alpha0;
                                num_C23=exp(-beta0*(R_s-1));
                                den_C21=beta0;
                                den_C22=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                den_C23=alpha0+R_s*beta0;

                                num_C31=temp_beta_1;
                                den_C31=beta0;
                                den_C32=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                num_C41=R_s-1;
                                num_C42=temp_beta_1;
                                den_C41=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                I_1_as=2*R_s*(num_A11/den_A11-num_A21/den_A21...
                                    -num_A31*num_A32/(den_A31*den_A32) +num_A41*num_A42/(den_A41*den_A42));
                                I_2_as=2*R_s*(num_B11/den_B11...
                                    -num_B21*num_B22*num_B23/(den_B21*den_B22*den_B23));
                                I_3_as=2*(num_C11*num_C12/(den_C11*den_C12)...
                                    +num_C21*num_C22*num_C23/(den_C21*den_C22*den_C23)...
                                    -num_C31/(den_C31*den_C32)+num_C41*num_C42/den_C41);

                                %                                 temp_asymp_out=temp_asymp_out-(-1)^i_m*(I_1_as+I_2_as+I_3_as)/(1/beta1(i_setrow));


                            end
                        end
                    end
                end   %if i_m==2


                if i_m==3
                    for i_1=1:no_relay-2
                        for i_2=i_1+1:no_relay-1
                            for i_3=i_2+1:no_relay
                                if i_1~=i_relay & i_2~=i_relay & i_3~=i_relay
                                    temp_beta_1=beta(chosen_set_append(i_1))+beta(chosen_set_append(i_2))...
                                        +beta(chosen_set_append(i_3));

                                    %----------------------------------------------------------------------
                                    %     beta_D=beta1+beta2;
                                    %     alpha_E=beta1+alpha1;

                                    B_1=beta_D(chosen_set_append(i_setrow,i_1))*beta0/(beta_D(chosen_set_append(i_setrow,i_1))-beta0);
                                    B_2=beta_D(chosen_set_append(i_setrow,i_1))*beta0/(beta0-beta_D(chosen_set_append(i_setrow,i_1)));
                                    num_A=B_1*exp(-beta0*(R_s-1));
                                    den_A1=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;
                                    den_A2=alpha0+R_s*beta0;

                                    num_B=B_2*exp(-beta_D(chosen_set_append(i_setrow,i_1))*(R_s-1));
                                    den_B1=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta_D(chosen_set_append(i_setrow,i_1));
                                    den_B2=alpha0+R_s*beta_D(chosen_set_append(i_setrow,i_1));

                                    num_C=num_A;
                                    den_C1=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;
                                    den_C2=alpha0+R_s*beta0;

                                    num_D=num_B;
                                    den_D1=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta_D(chosen_set_append(i_setrow,i_1));
                                    den_D2=alpha0+R_s*beta_D(chosen_set_append(i_setrow,i_1));

                                    A=num_A/(den_A1*den_A2);
                                    B=num_B/(den_B1*den_B2);
                                    C=num_C/(den_C1*den_C2);
                                    D=num_B/(den_D1*den_D2);


                                    E=B_1/(beta0*alpha0);
                                    F=B_2/(beta_D(i_relay)*alpha0);
                                    G=num_A/(beta0*den_A2);
                                    H=num_B/(beta_D(i_relay)*den_B2);

                                    I_1=R_s*alpha0*(A+B-C-D);
                                    I_2=R_s*temp_beta_1*alpha0*(C+D)/(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1);
                                    I_3=temp_beta_1*alpha0*(E+F-G-H)/(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1);

                                    temp_out_sel=temp_out_sel-(-1)^i_m*(I_1+I_2+I_3);
                                    %---------------------------------------------------

                                    num_A11=1;
                                    den_A11=alpha_E(chosen_set_append(i_setrow,i_1));

                                    num_A21=1;
                                    den_A21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                    num_A31=alpha0;
                                    num_A32=exp(-beta0*(R_s-1));
                                    den_A31=alpha0+R_s*beta0;
                                    den_A32=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;

                                    num_A41=alpha0;
                                    num_A42=exp(-beta0*(R_s-1));
                                    den_A41=alpha0+R_s*beta0;
                                    den_A42=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                    num_B11=alpha_E(chosen_set_append(i_setrow,i_1));
                                    den_B11=(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1).^2;

                                    num_B21=alpha0;
                                    num_B22=alpha_E(chosen_set_append(i_setrow,i_1));
                                    num_B23=exp(-beta0*(R_s-1));
                                    den_B21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                    den_B22=alpha0+R_s*beta0;
                                    den_B23=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                    num_C11=R_s;
                                    num_C12=temp_beta_1;
                                    den_C11=alpha0;
                                    den_C12=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                    num_C21=temp_beta_1;
                                    num_C22=alpha0;
                                    num_C23=exp(-beta0*(R_s-1));
                                    den_C21=beta0;
                                    den_C22=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                    den_C23=alpha0+R_s*beta0;

                                    num_C31=temp_beta_1;
                                    den_C31=beta0;
                                    den_C32=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                    num_C41=R_s-1;
                                    num_C42=temp_beta_1;
                                    den_C41=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                    I_1_as=2*R_s*(num_A11/den_A11-num_A21/den_A21...
                                        -num_A31*num_A32/(den_A31*den_A32) +num_A41*num_A42/(den_A41*den_A42));
                                    I_2_as=2*R_s*(num_B11/den_B11...
                                        -num_B21*num_B22*num_B23/(den_B21*den_B22*den_B23));
                                    I_3_as=2*(num_C11*num_C12/(den_C11*den_C12)...
                                        +num_C21*num_C22*num_C23/(den_C21*den_C22*den_C23)...
                                        -num_C31/(den_C31*den_C32)+num_C41*num_C42/den_C41);

                                    %                                     temp_asymp_out=temp_asymp_out-(-1)^i_m*(I_1_as+I_2_as+I_3_as)/(1/beta1(i_setrow));


                                end
                            end
                        end
                    end
                end   %if i_m==3

                if i_m==4
                    for i_1=1:no_relay-3
                        for i_2=i_1+1:no_relay-2
                            for i_3=i_2+1:no_relay-1
                                for i_4=i_3+1:no_relay
                                    if i_1~=i_relay & i_2~=i_relay & i_3~=i_relay & i_4~=i_relay
                                        temp_beta_1=beta(chosen_set_append(i_1))+beta(chosen_set_append(i_2))...
                                            +beta(chosen_set_append(i_3))+beta(chosen_set_append(i_4));

                                        %----------------------------------------------------------------------
                                        %     beta_D=beta1+beta2;
                                        %     alpha_E=beta1+alpha1;

                                        B_1=beta_D(chosen_set_append(i_setrow,i_1))*beta0/(beta_D(chosen_set_append(i_setrow,i_1))-beta0);
                                        B_2=beta_D(chosen_set_append(i_setrow,i_1))*beta0/(beta0-beta_D(chosen_set_append(i_setrow,i_1)));
                                        num_A=B_1*exp(-beta0*(R_s-1));
                                        den_A1=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;
                                        den_A2=alpha0+R_s*beta0;

                                        num_B=B_2*exp(-beta_D(chosen_set_append(i_setrow,i_1))*(R_s-1));
                                        den_B1=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta_D(chosen_set_append(i_setrow,i_1));
                                        den_B2=alpha0+R_s*beta_D(chosen_set_append(i_setrow,i_1));

                                        num_C=num_A;
                                        den_C1=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;
                                        den_C2=alpha0+R_s*beta0;

                                        num_D=num_B;
                                        den_D1=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta_D(chosen_set_append(i_setrow,i_1));
                                        den_D2=alpha0+R_s*beta_D(chosen_set_append(i_setrow,i_1));

                                        A=num_A/(den_A1*den_A2);
                                        B=num_B/(den_B1*den_B2);
                                        C=num_C/(den_C1*den_C2);
                                        D=num_B/(den_D1*den_D2);


                                        E=B_1/(beta0*alpha0);
                                        F=B_2/(beta_D(i_relay)*alpha0);
                                        G=num_A/(beta0*den_A2);
                                        H=num_B/(beta_D(i_relay)*den_B2);

                                        I_1=R_s*alpha0*(A+B-C-D);
                                        I_2=R_s*temp_beta_1*alpha0*(C+D)/(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1);
                                        I_3=temp_beta_1*alpha0*(E+F-G-H)/(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1);

                                        temp_out_sel=temp_out_sel-(-1)^i_m*(I_1+I_2+I_3);
                                        %---------------------------------------------------

                                        num_A11=1;
                                        den_A11=alpha_E(chosen_set_append(i_setrow,i_1));

                                        num_A21=1;
                                        den_A21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                        num_A31=alpha0;
                                        num_A32=exp(-beta0*(R_s-1));
                                        den_A31=alpha0+R_s*beta0;
                                        den_A32=alpha_E(chosen_set_append(i_setrow,i_1))+R_s*beta0;

                                        num_A41=alpha0;
                                        num_A42=exp(-beta0*(R_s-1));
                                        den_A41=alpha0+R_s*beta0;
                                        den_A42=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                        num_B11=alpha_E(chosen_set_append(i_setrow,i_1));
                                        den_B11=(alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1).^2;

                                        num_B21=alpha0;
                                        num_B22=alpha_E(chosen_set_append(i_setrow,i_1));
                                        num_B23=exp(-beta0*(R_s-1));
                                        den_B21=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                        den_B22=alpha0+R_s*beta0;
                                        den_B23=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1+R_s*beta0;

                                        num_C11=R_s;
                                        num_C12=temp_beta_1;
                                        den_C11=alpha0;
                                        den_C12=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                        num_C21=temp_beta_1;
                                        num_C22=alpha0;
                                        num_C23=exp(-beta0*(R_s-1));
                                        den_C21=beta0;
                                        den_C22=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;
                                        den_C23=alpha0+R_s*beta0;

                                        num_C31=temp_beta_1;
                                        den_C31=beta0;
                                        den_C32=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                        num_C41=R_s-1;
                                        num_C42=temp_beta_1;
                                        den_C41=alpha_E(chosen_set_append(i_setrow,i_1))+temp_beta_1;

                                        I_1_as=2*R_s*(num_A11/den_A11-num_A21/den_A21...
                                            -num_A31*num_A32/(den_A31*den_A32) +num_A41*num_A42/(den_A41*den_A42));
                                        I_2_as=2*R_s*(num_B11/den_B11...
                                            -num_B21*num_B22*num_B23/(den_B21*den_B22*den_B23));
                                        I_3_as=2*(num_C11*num_C12/(den_C11*den_C12)...
                                            +num_C21*num_C22*num_C23/(den_C21*den_C22*den_C23)...
                                            -num_C31/(den_C31*den_C32)+num_C41*num_C42/den_C41);

                                        %                                         temp_asymp_out=temp_asymp_out-(-1)^i_m*(I_1_as+I_2_as+I_3_as)/(1/beta1(i_setrow));

                                    end
                                end

                            end
                        end
                    end
                end   %if i_m==4
            end
            %           end
            %                 end  %i_relay=1:no_relay

            %----------------------------------------------------
            %
            %         end    % i_m
            %         %         temp_sum_bino=0;
            %         %         for i_bino=0:1:no_relay
            %         %             temp_sum_bino=temp_sum_bino+nchoosek(no_relay,i_bino) *(R_s-1)^i_bino*R_s^(no_relay-i_bino)*gamma(no_relay-i_bino+1)...
            %         %                 /(alpha(i_relay)^(no_relay-i_bino));
            %         %
            %         end
            %         %
            %         %         temp_asymp_out=temp_asymp_out+ prod(beta) * temp_sum_bino/no_relay;
            %

            %         end
            %----------------------------------------------------


            %     sec_outg_asymp(i_snr)=2*beta1*((R_s-1)+alpha0*alpha1*exp(-beta0*(R_s-1))...
            %              /(beta0*(alpha0+R_s*beta0)*(alpha1+R_s*beta0)));

            %     sec_outg(i_snr)=temp_out_sel;
            %     sec_outg_asymp(i_snr)= temp_asymp_out
            %     i_snr=i_snr+1;
            prod_sel=1;
            for i_setsel=1:i_relay
                prod_sel=prod_sel*exp(-beta1(chosen_set(i_setrow,i_setsel))*gamma_th);
            end
            comp_set=setdiff(relay_set,chosen_set(i_setrow,:));
            prod_notsel=1;
            for i_comset=1:length(comp_set)
                prod_notsel=prod_notsel*(1-exp(-beta1(comp_set(i_comset))*gamma_th));
            end
            temp_out=temp_out+prod_sel*prod_notsel*temp_out_sel;
            %             end
        end  %i_setrow
        final_outg=final_outg+temp_out;
        prob_nosel=prob_nosel*(1-exp(-beta1(i_relay)*gamma_th));
    end   %i_relay
    
    sec_outg(i_snr)=final_outg+prob_nosel*(1-alpha0*exp(-(R_s-1)*beta0)/(alpha0+R_s*beta0));
%     end
    
    i_snr=i_snr+1;
end   %i_snr
% =========================================================
semilogy(snr_dB,sec_outg,'-or');
hold on;
% semilogy(snr_dB,sec_outg_asymp,'-k');
% hold on;
%  --------------------------------------------------------------
%  save('ANA_RSEL_MAXE_N_1_R_1_0_E1_3.mat', 'snr_dB', 'sec_outg','sec_outg_asymp','no_relay');
% save('ANA_RSEL_MAXE_N_2_R_1_0_E1_0_3.mat', 'snr_dB', 'sec_outg','no_relay');
%  save('ANA_RSEL_MAXE_N_3_R_1_0_E1_0_3_6mat', 'snr_dB', 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_4_R_1_0_E1_0_3_6_9.mat','snr_dB','sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_4_R_1_0_E1_3_6_9_12.mat','snr_dB','sec_outg','no_relay');

% save('ANA_RSEL_MAXE_N_4_R_0_E1_0_3_6_9.mat','snr_dB','sec_outg','no_relay');

% save('ANA_RSEL_MAXE_N_2_R_1_0_E1_3_3.mat', 'snr_dB', 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_4_R_1_0_E1_3_3_3_3.mat', 'snr_dB', 'sec_outg','no_relay');

% save('ANA_RSEL_MAXE_N_2_R_1_0_E1_0_3_B.mat', 'snr_dB', 'sec_outg','no_relay');
%  save('ANA_RSEL_MAXE_N_3_R_1_0_E1_0_3_6_B.mat', 'snr_dB', 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_4_R_1_0_E1_0_3_6_9_B.mat','snr_dB','sec_outg','no_relay');

% save('ANA_RSEL_MAXE_N_1_R_1_0_E1_3_B.mat', 'snr_dB',...
% 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_2_R_1_0_E1_3_3_B.mat', 'snr_dB', 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_3_R_1_0_E1_3_3_3_B.mat', 'snr_dB',...
% 'sec_outg','no_relay');
% save('ANA_RSEL_MAXE_N_4_R_1_0_E1_3_3_3_3_B.mat', 'snr_dB',...
% 'sec_outg','no_relay');

%  --------------------------------------------------------------

% grid on;
xlabel('Total SNR');
ylabel('P_{o}(R_s)');
% axis([min(snr_dB)  max(snr_dB)  10^-4 1]);
% axis([min(snr_dB)  40  10^-4 1]);
% axis([min(snr_dB)  40  10^-3 1]);

% print -deps ANA_RSEL_TRAD_N_4_R_1_1_E_3
% =========================================================
