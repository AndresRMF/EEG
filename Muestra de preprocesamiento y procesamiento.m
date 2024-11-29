signals=[];
eegdataa= load('eegData3.mat');
if app.ShirtsignalsCheckBox.Value
    signals = {app.ab1,app.ab3,app.ab5,app.an1,app.an3,app.an5,app.gt1,app.gt3,app.gt5,...
               app.g1, app.g3, app.g5, app.gy1,app.gy3,app.gy5,app.kb1,app.kb3,app.kb5,app.mh1,...
                       app.mh3,app.mh5,app.mt1,app.mt3,app.mt5,app.ps1,app.ps3,app.ps5,app.pd1,...
                       app.pd3,app.pd5,app.re1,app.re3,app.re5,app.rj1,app.rj3,app.rj5,app.rb1,...
                       app.rb3,app.rb5,app.ra1,app.ra3,app.ra5,app.rs1,app.rs3,app.rs5,app.rk1,...
                       app.rk3,app.rk5,app.sc1,app.sc3,app.sc5,app.sd1,app.sd3,app.sd5,app.sm1,...
                       app.sm3,app.sm5,app.ss1,app.ss3,app.ss5,app.tq1,app.tq3,app.tq5,app.vp1,...
                       app.vp3,app.vp5,app.vj1,app.vj3,app.vj5,app.vi1,app.vi3,app.vi5,app.vr1,...
                       app.vr3,app.vr5};
elseif app.SweatersignalsCheckBox.Value
    signals = {app.ab2,app.ab4,app.ab6,app.an2,app.an4,app.an6,app.gt2,app.gt4,app.gt6,...
               app.g2, app.g4, app.g6, app.gy2,app.gy4,app.gy6,app.kb2,app.kb4,app.kb6,...
               app.mh2,app.mh4,app.mh6,app.mt2,app.mt4,app.mt6,app.ps2,app.ps4,app.ps6,...
               app.pd2,app.pd4,app.pd6,app.re2,app.re4,app.re6,app.rj2,app.rj4,app.rj6,...
               app.rb2,app.rb4,app.rb6,app.ra2,app.ra4,app.ra6,app.rs2,app.rs4,app.rs6,...
               app.rk2,app.rk4,app.rk6,app.sc2,app.sc4,app.sc6,app.sd2,app.sd4,app.sd6,...
               app.sm2,app.sm4,app.sm6,app.ss2,app.ss4,app.ss6,app.tq2,app.tq4,app.tq6,...
               app.vp2,app.vp4,app.vp6,app.vj2,app.vj4,app.vj6,app.vi2,app.vi4,app.vi6,...
               app.vr2,app.vr4,app.vr6};
elseif app.BothsignalsCheckBox.Value
    signals = {app.ab1, app.ab2, app.ab3, app.ab4, app.ab5,app.ab6,app.an1,app.an2,app.an3,app.an4,app.an5,app.an6,app.gt1,app.gt2,app.gt3,app.gt4,app.gt5,app.gt6,...
       app.g1,app.g2,app.g3,app.g4,app.g5,app.g6,app.gy1,app.gy2,app.gy3,app.gy4,app.gy5,app.gy6,app.kb1,app.kb2,app.kb3,app.kb4,app.kb5,app.kb6,app.mh1,...
       app.mh2,app.mh3,app.mh4,app.mh5,app.mh6,app.mt1,app.mt2,app.mt3,app.mt4,app.mt5,app.mt6,app.ps1,app.ps2,app.ps3,app.ps4,app.ps5,app.ps6,app.pd1,...
       app.pd2,app.pd3,app.pd4,app.pd5,app.pd6,app.re1,app.re2,app.re3,app.re4,app.re5,app.re6,app.rj1,app.rj2,app.rj3,app.rj4,app.rj5,app.rj6,app.rb1,...
       app.rb2,app.rb3,app.rb4,app.rb5,app.rb6,app.ra1,app.ra2,app.ra3,app.ra4,app.ra5,app.ra6,app.rs1,app.rs2,app.rs3,app.rs4,app.rs5,app.rs6,app.rk1,...
       app.rk2,app.rk3,app.rk4,app.rk5,app.rk6,app.sc1,app.sc2,app.sc3,app.sc4,app.sc5,app.sc6,app.sd1,app.sd2,app.sd3,app.sd4,app.sd5,app.sd6,app.sm1,...
       app.sm2,app.sm3,app.sm4,app.sm5,app.sm6,app.ss1,app.ss2,app.ss3,app.ss4,app.ss5,app.ss6,app.tq1,app.tq2,app.tq3,app.tq4,app.tq5,app.tq6,app.vp1,...
       app.vp2,app.vp3,app.vp4,app.vp5,app.vp6,app.vj1,app.vj2,app.vj3,app.vj4,app.vj5,app.vj6,app.vi1,app.vi2,app.vi3,app.vi4,app.vi5,app.vi6,app.vr1,...
       app.vr2,app.vr3,app.vr4,app.vr5,app.vr6};
elseif app.AlldatasetCheckBox.Value
    signals = eegdataa.eeg_Data;
end

numSignal = length(signals);
matrizCaracteristicas = []; % Inicializa la matriz de características
%eegSeleccionada = [];
muestrasPorVentana = app.WindowLengthEditField_2.Value;  % Tamaño de cada ventana en muestras
solapamientoMuestras = app.OverlappingEditField.Value; % Muestras de solapamiento entre ventanas

totalVentanasEstimadas = sum(arrayfun(@(x) max(0, floor((size(signals{x}, 1) - muestrasPorVentana) / (muestrasPorVentana - solapamientoMuestras) + 1)), 1:numSignal));
indiceMatriz = 1; % Índice para llenar la matriz de características

% Crea la waitbar
h = waitbar(0, 'Generating the feature matrix...');

for idx = 1:numSignal

    eegSeleccionada = signals{idx};

    [numMuestras, numCanales] = size(eegSeleccionada);

    % Aplicar el filtro Pasa Banda si está seleccionado
    if app.UseCheckBox_2.Value
        f_low = app.LowcutofffrequencyEditField.Value; % Frecuencia de corte baja
        f_high = app.HighcutofffrequencyEditField.Value; % Frecuencia de corte alta
        eegSeleccionada = filtrarPasaBanda(app,eegSeleccionada, app.Fs, f_low, f_high);

    end

    % Aplicar el filtro Notch si está seleccionado
    if app.UseCheckBox.Value
        f0 = app.FrequencytodeleteEditField.Value; % Frecuencia que se desea eliminar
        bw = 1; % Ancho de banda del filtro Notch

        eegSeleccionada = aplicarFiltroNotch(app,eegSeleccionada, app.Fs, f0, bw);
    end

    % Aplicar el filtro Savitzky-Golay si está seleccionado
    if app.UseCheckBox_3.Value
        orden = app.PolynomialOrderEditField.Value;
        sizef = app.WindowLengthEditField.Value;

        eegSeleccionada = aplicarSavitzkyGolay(app,eegSeleccionada,orden, sizef);

    end

    if app.IntervalsCheckBox.Value

        % Bucle para extraer ventanas con solapamiento
        for indiceInicio = 1:(muestrasPorVentana - solapamientoMuestras):numMuestras
            indiceFin = indiceInicio + muestrasPorVentana - 1;

            % Verifica si el índice final excede el número de muestras disponibles
            if indiceFin > numMuestras
                break;  % Sale del loop si no hay suficientes muestras para completar la ventana
            end

            % Extraer ventana de la señal
            eegVentana = eegSeleccionada(indiceInicio:indiceFin, :);

            % % Crear y aplicar ventana de Hamming
            hanningWindow = hann(length(eegVentana));
            eegVentana = eegVentana .* hanningWindow;

            caracteristicasActual=[];

            fs = app.Fs;

            % Asignación de etiquetas a cada columna de la señal EEG
            etiquetasCanales = {'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};

            % Definición de los grupos de canales con etiquetas
            gruposCanales = containers.Map();
            gruposCanales('Frontal') = {'FC6', 'F4', 'F8', 'AF4', 'AF3', 'F7', 'F3', 'FC5'};
            gruposCanales('Temporal') = {'T7', 'T8'};
            gruposCanales('Occipital') = {'O1', 'O2'};
            gruposCanales('Parietal') = {'P7', 'P8'};

            % Verifica cada grupo y agrega los índices de los canales seleccionados
            indicesCanalesSeleccionados = [];
            if app.FrontalChannelsCheckBox.Value
                indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Frontal'), etiquetasCanales)];
            end
            if app.TemporalChannelsCheckBox.Value
                indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Temporal'), etiquetasCanales)];
            end
            if app.OccipitalChannelsCheckBox.Value
                indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Occipital'), etiquetasCanales)];
            end
            if app.ParientalChannelsCheckBox.Value
                indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Parietal'), etiquetasCanales)];
            end

            % Ordena los índices para mantener el orden original de los canales
            indicesCanalesSeleccionados = sort(indicesCanalesSeleccionados);

            % Extrae los canales seleccionados de la matriz EEG
            eegVentana = eegVentana(:, indicesCanalesSeleccionados);

            [numSamples, numChannels] = size(eegVentana);

            if app.FrequencyFeaturesCheckBox.Value

                % Inicializa las variables de salida
                delta = 0;
                theta = 0;
                alpha = 0;
                beta = 0;
                gamma = 0;

                % Inicialización del vector de características global
                globalFeaturesVector = [];

                for ch = 1:numChannels
                    data = eegVentana(:, ch);
                    fftData = fft(data);
                    P2 = abs(fftData / numSamples);
                    P1 = P2(1:numSamples / 2 + 1)';
                    P1(2:end) = 2 * P1(2:end);
                    f = fs * (0:(numSamples / 2)) / numSamples;

                    channelFeatures = [];

                    % Calcula el band power para las bandas definidas
                    delta = delta + sum(P1(f >= 0.5 & f <= 4));
                    theta = theta + sum(P1(f >= 4 & f <= 8));
                    alpha = alpha + sum(P1(f >= 8 & f <= 13));
                    beta = beta + sum(P1(f >= 13 & f <= 30));
                    gamma = gamma + sum(P1(f >= 30 & f <= 50));

                    % Verificar cada checkbox y recopilar las bandas correspondientes
                    if app.Gamma30hzCheckBox.Value
                        channelFeatures = [channelFeatures, gamma];
                    end
                    if app.Beta1330hzCheckBox.Value
                        channelFeatures = [channelFeatures, beta];
                    end
                    if app.Alpha813hzCheckBox.Value
                        channelFeatures = [channelFeatures, alpha];
                    end
                    if app.Theta48hzCheckBox.Value
                        channelFeatures = [channelFeatures, theta];
                    end
                    if app.Delta054hzCheckBox.Value
                        channelFeatures = [channelFeatures, delta];
                    end

                    spectralEntropy = -sum((P1 / sum(P1)) .* log2((P1 / sum(P1)) + eps));

                    % Agregar solo las características calculadas por canal
                    channelFeatures = [channelFeatures, spectralEntropy];

                    globalFeaturesVector = [globalFeaturesVector, channelFeatures];
                end
                caracteristicasActual = [caracteristicasActual, globalFeaturesVector];
            end

            if app.TimeFreqFeaturesCheckBox.Value
                % Inicialización del vector de características global
                globalFeaturesVector = [];

                % Bucle para cada canal
                for ch = 1:numChannels
                    data = eegVentana(:, ch);

                    % Transformada de Wavelet
                    [C, L] = wavedec(data, 5, 'db4');
                    cA = appcoef(C, L, 'db4');  % Coeficientes de aproximación al nivel más bajo
                    cD1 = detcoef(C, L, 1);     % Coeficientes de detalle nivel 1 (Gamma alta)
                    cD2 = detcoef(C, L, 2);     % Coeficientes de detalle nivel 2 (Beta-Gamma baja)
                    cD3 = detcoef(C, L, 3);     % Coeficientes de detalle nivel 3 (Alpha)
                    cD4 = detcoef(C, L, 4);     % Coeficientes de detalle nivel 4 (Theta)
                    cD5 = detcoef(C, L, 5);     % Coeficientes de detalle nivel 5 (Delta baja)

                    % Calculando la energía para cada banda de interés
                    energyDelta = sum(cD5.^2);  % Considerar cD5 como Delta
                    energyTheta = sum(cD4.^2);  % Considerar cD4 como Theta
                    energyAlpha = sum(cD3.^2);  % Considerar cD3 como Alpha
                    energyBeta = sum(cD2.^2);   % Considerar cD2 como Beta
                    energyGamma = sum(cD1.^2);  % Considerar cD1 como Gamma

                    selectedEnergy = [];
                    % Verificar cada checkbox y recopilar las bandas correspondientes
                    if app.Gamma30hzCheckBox.Value
                        selectedEnergy = [selectedEnergy, energyGamma];
                    end
                    if app.Beta1330hzCheckBox.Value
                        selectedEnergy = [selectedEnergy, energyBeta];
                    end
                    if app.Alpha813hzCheckBox.Value
                        selectedEnergy = [selectedEnergy, energyAlpha];
                    end
                    if app.Theta48hzCheckBox.Value
                        selectedEnergy = [selectedEnergy, energyTheta];
                    end
                    if app.Delta054hzCheckBox.Value
                        selectedEnergy = [selectedEnergy, energyDelta];
                    end

                    % Calculando características Hjorth
                    activity = var(data);
                    mobility = sqrt(var(diff(data)) / activity);

                    % Construir vector de características para el canal actual
                    channelFeatures = [activity, mobility, selectedEnergy];

                    % Concatenar al vector de características global
                    globalFeaturesVector = [globalFeaturesVector, channelFeatures];
                end
                caracteristicasActual=[caracteristicasActual, globalFeaturesVector];
            end

            if app.AlphaFeaturesCheckBox.Value
                featuresa=alphaf(app,eegVentana);
                caracteristicasActual=[caracteristicasActual,featuresa];
            end

            if app.IntervalsCheckBox_4.Value

                if app.FFTCheckBox.Value
                    if app.MeanCheckBox.Value
                        [promedioDelta, promedioTheta, promedioAlpha, promedioBeta, promedioGamma] = calcularPromediosPSD(app,eegVentana);
                    else
                        [promedioDelta, promedioTheta, promedioAlpha, promedioBeta, promedioGamma] = PSD(app,eegVentana);
                    end

                elseif app.WVCheckBox.Value
                    [promedioGamma, promedioBeta, promedioAlpha, promedioTheta, promedioDelta] = discwav(app,eegSeleccionada);
                    promedioGamma = reshape(promedioGamma, 1, []);
                    promedioBeta=reshape(promedioBeta, 1, []);
                    promedioAlpha=reshape(promedioAlpha, 1, []);
                    promedioTheta=reshape(promedioTheta, 1, []);
                    promedioDelta=reshape(promedioDelta, 1, []);

                    if app.MeanCheckBox.Value
                        promedioGamma = mean(promedioGamma);
                        promedioBeta = mean(promedioBeta);
                        promedioAlpha = mean(promedioAlpha);
                        promedioTheta = mean(promedioTheta);
                        promedioDelta = mean(promedioDelta);
                    end
                end
                
                bandasSeleccionadas = []; % Vector para almacenar las bandas seleccionadas

                % Verificar cada checkbox y recopilar las bandas correspondientes
                if app.Gamma30hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioGamma];
                end
                if app.Beta1330hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioBeta];
                end
                if app.Alpha813hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioAlpha];
                end
                if app.Theta48hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioTheta];
                end
                if app.Delta054hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioDelta];
                end
                caracteristicasActual=[caracteristicasActual,bandasSeleccionadas];
            end

            % Añade las características de esta ventana a la matriz general
            matrizCaracteristicas(indiceMatriz, :) = caracteristicasActual;
            indiceMatriz = indiceMatriz + 1; % Incrementa el índice para la próxima ventana

             % Actualizar la waitbar cada 5 iteraciones
            if mod(idx, 10) == 0
                waitbar(idx / numSignal, h, sprintf('Processing: %d %%', floor(idx / numSignal * 100)));
            end
        end
    else
        caracteristicasActual=[];

        fs = app.Fs;

        % Definir los índices de los canales F3 y F4 en el orden proporcionado
        indice_F3 = 3;  % Canal 'F3' es el tercer canal
        indice_F4 = 12; % Canal 'F4' es el duodécimo canal

        % Extraer los datos de los canales F3 y F4
        datos_F3 = eegSeleccionada(:, indice_F3);
        datos_F4 = eegSeleccionada(:, indice_F4);

        % Asignación de etiquetas a cada columna de la señal EEG
        etiquetasCanales = {'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};

        % Definición de los grupos de canales con etiquetas
        gruposCanales = containers.Map();
        gruposCanales('Frontal') = {'FC6', 'F4', 'F8', 'AF4', 'AF3', 'F7', 'F3', 'FC5'};
        gruposCanales('Temporal') = {'T7', 'T8'};
        gruposCanales('Occipital') = {'O1', 'O2'};
        gruposCanales('Parietal') = {'P7', 'P8'};

        % Verifica cada grupo y agrega los índices de los canales seleccionados
        indicesCanalesSeleccionados = [];
        if app.FrontalChannelsCheckBox.Value
            indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Frontal'), etiquetasCanales)];
        end
        if app.TemporalChannelsCheckBox.Value
            indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Temporal'), etiquetasCanales)];
        end
        if app.OccipitalChannelsCheckBox.Value
            indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Occipital'), etiquetasCanales)];
        end
        if app.ParientalChannelsCheckBox.Value
            indicesCanalesSeleccionados = [indicesCanalesSeleccionados, encontrarIndicesCanales(app,gruposCanales('Parietal'), etiquetasCanales)];
        end

        % Ordena los índices para mantener el orden original de los canales
        indicesCanalesSeleccionados = sort(indicesCanalesSeleccionados);

        % Extrae los canales seleccionados de la matriz EEG
        eegSeleccionada = eegSeleccionada(:, indicesCanalesSeleccionados);

        [numSamples, numChannels] = size(eegSeleccionada);

        if app.FrequencyFeaturesCheckBox.Value
            % Definir las bandas de frecuencia
            %bandas = [0.5 4; 4 8; 8 13; 13 30; 30 50];  % Delta, Theta, Alpha, Beta, Gamma
            %numBands = size(bandas, 1);

            % Inicializa las variables de salida
            delta = 0;
            theta = 0;
            alpha = 0;
            beta = 0;
            gamma = 0;

            % Inicialización del vector de características global
            globalFeaturesVector = [];

            for ch = 1:numChannels
                data = eegSeleccionada(:, ch);
                fftData = fft(data);
                P2 = abs(fftData / numSamples);
                P1 = P2(1:numSamples / 2 + 1)';
                P1(2:end) = 2 * P1(2:end);
                f = fs * (0:(numSamples / 2)) / numSamples;

                channelFeatures = [];

                % for b = 1:numBands
                %     freqRange = bandas(b, :);
                %     bandIndices = f >= freqRange(1) & f <= freqRange(2);
                %     bandPower = sum(P1(bandIndices));
                %     channelFeatures = [channelFeatures, bandPower];
                % end

                % Calcula el band power para las bandas definidas
                delta = delta + sum(P1(f >= 0.5 & f <= 4));
                theta = theta + sum(P1(f >= 4 & f <= 8));
                alpha = alpha + sum(P1(f >= 8 & f <= 13));
                beta = beta + sum(P1(f >= 13 & f <= 30));
                gamma = gamma + sum(P1(f >= 30 & f <= 50));

                % Verificar cada checkbox y recopilar las bandas correspondientes
                if app.Gamma30hzCheckBox.Value
                    channelFeatures = [channelFeatures, gamma];
                end
                if app.Beta1330hzCheckBox.Value
                    channelFeatures = [channelFeatures, beta];
                end
                if app.Alpha813hzCheckBox.Value
                    channelFeatures = [channelFeatures, alpha];
                end
                if app.Theta48hzCheckBox.Value
                    channelFeatures = [channelFeatures, theta];
                end
                if app.Delta054hzCheckBox.Value
                    channelFeatures = [channelFeatures, delta];
                end
                
                spectralEntropy = -sum((P1 / sum(P1)) .* log2((P1 / sum(P1)) + eps));

                %spectralCentroid = sum(f .* P1) / sum(P1);

                %spectralSpread = sqrt(sum((f - spectralCentroid).^2 .* P1) / sum(P1));

                % Agregar solo las características calculadas por canal
                channelFeatures = [channelFeatures, spectralEntropy];

                globalFeaturesVector = [globalFeaturesVector, channelFeatures];
            end
            caracteristicasActual = [caracteristicasActual, globalFeaturesVector];
        end

        if app.TimeFreqFeaturesCheckBox.Value
            % Inicialización del vector de características global
            globalFeaturesVector = [];

            % Bucle para cada canal
            for ch = 1:numChannels
                data = eegSeleccionada(:, ch);
                % Transformada de Wavelet
                [C, L] = wavedec(data, 5, 'db4');
                cA = appcoef(C, L, 'db4');  % Coeficientes de aproximación al nivel más bajo
                cD1 = detcoef(C, L, 1);     % Coeficientes de detalle nivel 1 (Gamma alta)
                cD2 = detcoef(C, L, 2);     % Coeficientes de detalle nivel 2 (Beta-Gamma baja)
                cD3 = detcoef(C, L, 3);     % Coeficientes de detalle nivel 3 (Alpha)
                cD4 = detcoef(C, L, 4);     % Coeficientes de detalle nivel 4 (Theta)
                cD5 = detcoef(C, L, 5);     % Coeficientes de detalle nivel 5 (Delta baja)

                % Calculando la energía para cada banda de interés
                energyDelta = sum(cD5.^2);  % Considerar cD5 como Delta
                energyTheta = sum(cD4.^2);  % Considerar cD4 como Theta
                energyAlpha = sum(cD3.^2);  % Considerar cD3 como Alpha
                energyBeta = sum(cD2.^2);   % Considerar cD2 como Beta
                energyGamma = sum(cD1.^2);  % Considerar cD1 como Gamma

                selectedEnergy = [];
                 % Verificar cada checkbox y recopilar las bandas correspondientes
                if app.Gamma30hzCheckBox.Value
                    selectedEnergy = [selectedEnergy, energyGamma];
                end
                if app.Beta1330hzCheckBox.Value
                    selectedEnergy = [selectedEnergy, energyBeta];
                end
                if app.Alpha813hzCheckBox.Value
                    selectedEnergy = [selectedEnergy, energyAlpha];
                end
                if app.Theta48hzCheckBox.Value
                    selectedEnergy = [selectedEnergy, energyTheta];
                end
                if app.Delta054hzCheckBox.Value
                    selectedEnergy = [selectedEnergy, energyDelta];
                end

                % Calculando características Hjorth
                activity = var(data);
                mobility = sqrt(var(diff(data)) / activity);

                % Construir vector de características para el canal actual
                channelFeatures = [activity, mobility, selectedEnergy];

                % Concatenar al vector de características global
                globalFeaturesVector = [globalFeaturesVector, channelFeatures];
            end
            caracteristicasActual=[caracteristicasActual, globalFeaturesVector];
        end
        %disp(['Tamaño de caracteristicasActual después de la time-freq: ', num2str(length(caracteristicasActual))]);

        if app.AWCheckBox.Value
            aw_number=AW(app, datos_F3, datos_F4);
            caracteristicasActual=[caracteristicasActual, aw_number];
        end
        %disp(['Tamaño de caracteristicasActual después de la aw: ', num2str(length(caracteristicasActual))]);

        if app.AlphaFeaturesCheckBox.Value
            featuresa=alphaf(app,eegSeleccionada);
            caracteristicasActual=[caracteristicasActual,featuresa];
        end
        %disp(['Tamaño de caracteristicasActual después de la af: ', num2str(length(caracteristicasActual))]);

        if app.IntervalsCheckBox_4.Value

            if app.FFTCheckBox.Value
                if app.MeanCheckBox.Value
                    [promedioDelta, promedioTheta, promedioAlpha, promedioBeta, promedioGamma] = calcularPromediosPSD(app,eegSeleccionada);
                else
                    [promedioDelta, promedioTheta, promedioAlpha, promedioBeta, promedioGamma] = PSD(app,eegSeleccionada);
                end

                bandasSeleccionadas = []; % Vector para almacenar las bandas seleccionadas

                % Verificar cada checkbox y recopilar las bandas correspondientes
                if app.Gamma30hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioGamma];
                end
                if app.Beta1330hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioBeta];
                end
                if app.Alpha813hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioAlpha];
                end
                if app.Theta48hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioTheta];
                end
                if app.Delta054hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioDelta];
                end
                caracteristicasActual=[caracteristicasActual,bandasSeleccionadas];

            elseif app.WVCheckBox.Value
                if app.MeanCheckBox.Value
                    [promedioGamma, promedioBeta, promedioAlpha, promedioTheta, promedioDelta] = promdiscwav(app,eegSeleccionada);
                else
                    [promedioGamma, promedioBeta, promedioAlpha, promedioTheta, promedioDelta] = discwav(app,eegSeleccionada);
                    promedioGamma = reshape(promedioGamma, 1, []);
                    promedioBeta=reshape(promedioBeta, 1, []);
                    promedioAlpha=reshape(promedioAlpha, 1, []);
                    promedioTheta=reshape(promedioTheta, 1, []);
                    promedioDelta=reshape(promedioDelta, 1, []);
                end
                bandasSeleccionadas = []; % Vector para almacenar las bandas seleccionadas

                % Verificar cada checkbox y recopilar las bandas correspondientes
                if app.Gamma30hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioGamma];
                end
                if app.Beta1330hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioBeta];
                end
                if app.Alpha813hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioAlpha];
                end
                if app.Theta48hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioTheta];
                end
                if app.Delta054hzCheckBox.Value
                    bandasSeleccionadas = [bandasSeleccionadas, promedioDelta];
                end
                caracteristicasActual=[caracteristicasActual,bandasSeleccionadas];
            end
            %disp(['Tamaño de caracteristicasActual después de fft: ', num2str(length(caracteristicasActual))]);


        end
        %disp(['Tamaño de caracteristicasActual después de ejecutar el codigo: ', num2str(length(caracteristicasActual))]);

        % Añade el vector de características a la matriz
        matrizCaracteristicas(idx, :) = caracteristicasActual;

        % Actualizar la waitbar cada 5 iteraciones
        if mod(idx, 10) == 0
            waitbar(idx / numSignal, h, sprintf('Processing: %d %%', floor(idx / numSignal * 100)));
        end

    end
    
end

% Cierra la waitbar después de completar el bucle
waitbar(1, h, 'Feature extraction complete.');
close(h); % Cierra la waitbar

if app.NormalizedMinMaxCheckBox.Value

    % Obtener el nombre de la matriz del campo de texto
    matrixName = app.MatrixNameEditField.Value;
    % Verificar si el nombre de la matriz está vacío y asignar un nombre por defecto si es necesario
    if isempty(matrixName)
        matrixName = 'EEGData'; % Nombre por defecto
    end
    % Reemplazar caracteres no válidos en el nombre del archivo
    validMatrixName = matlab.lang.makeValidName(matrixName);

    % Solicitar al usuario que seleccione la carpeta donde guardar el archivo
    folder = '/MATLAB Drive/Tesis -Neuromarketing/Matrix';

    % Definir el nombre del archivo
    filename = fullfile(folder, [validMatrixName, '.mat']);

    % Guardar los datos EEG en un archivo MAT
    save(filename, 'EEG_data_3D'); % Reemplaza 'EEG_data_3D' con el nombre de tu variable de datos EEG 3D


else
    % Obtener el nombre de la matriz del campo de texto
    matrixName = app.MatrixNameEditField.Value;
    % Verificar si el nombre de la matriz está vacío y asignar un nombre por defecto si es necesario
    if isempty(matrixName)
        matrixName = 'EEGFeatures'; % Nombre por defecto
    end
    % Reemplazar caracteres no válidos en el nombre del archivo
    validMatrixName = matlab.lang.makeValidName(matrixName);

    % Solicitar al usuario que seleccione la carpeta donde guardar el archivo
    folder = '/MATLAB Drive/Tesis -Neuromarketing/Matrix';

    % Definir el nombre del archivo
    filename = fullfile(folder, [validMatrixName, '.csv']);

    % Guardar las características en un archivo CSV
    writematrix(matrizCaracteristicas, filename);

end