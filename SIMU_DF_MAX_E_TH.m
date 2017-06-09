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

    % ---------------creating channel------------------
    channel_0 = sqrt(1/beta0).*abs(sqrt(0.5)*(randn(1,no_simu)+i*randn(1,no_simu)));
    channel_e_0 = sqrt(1/alpha0).*abs(sqrt(0.5)*(randn(1,no_simu)+i*randn(1,no_simu)));

    % % % % % %   -----------------------------------------
    for i_relay=1:no_relay
        channel_1(i_relay,:) = sqrt(1/beta1(i_relay)).*abs(sqrt(0.5)*(randn(1,no_simu)+i*randn(1,no_simu)));
        channel_2(i_relay,:) = sqrt(1/beta2(i_relay)).*abs(sqrt(0.5)*(randn(1,no_simu)+i*randn(1,no_simu)));
        channel_e_1(i_relay,:) = sqrt(1/alpha1(i_relay)).*abs(sqrt(0.5)*(randn(1,no_simu)+i*randn(1,no_simu)));
    end

    no_df_direct=0;

    %%%%%--------------checking condition------------------------
    for i_simu=1:no_simu
        for i_relay=1:no_relay
            if(channel_1(i_relay,i_simu)^2>gamma_th)
                end_SNR_main(i_relay) = channel_0(i_simu)^2+channel_2(i_relay,i_simu)^2;
                end_SNR_evsd(i_relay) = channel_e_0(i_simu)^2+channel_e_1(i_relay,i_simu)^2;
            else
                end_SNR_main(i_relay)=channel_0(i_simu)^2;
                end_SNR_evsd(i_relay)=channel_e_0(i_simu)^2;
            end
        end

        %select max from R-E relays
        [max_capa,index_max] = max(channel_e_1(:,i_simu).^2,[],1);

        if (1+end_SNR_main(index_max))/(1+end_SNR_evsd(index_max))<R_s
            no_df_direct=no_df_direct+1;
        end
    end
    % end

    sec_outg(i_snr)=no_df_direct/no_simu;
    i_snr=i_snr+1;
end
% =========================================================
semilogy(snr_dB,sec_outg,'-+k');
hold on;
xlabel('Total SNR');
ylabel('P_{o}(R_s)');
axis([min(snr_dB)  max(snr_dB)  10^-4 1]);
