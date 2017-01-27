function varargout = SimApp(varargin)
% SIMAPP MATLAB code for SimApp.fig
%      SIMAPP, by itself, creates a new SIMAPP or raises the existing
%      singleton*.
%
%      H = SIMAPP returns the handle to a new SIMAPP or the handle to
%      the existing singleton*.
%
%      SIMAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMAPP.M with the given input arguments.
%
%      SIMAPP('Property','Value',...) creates a new SIMAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimApp

% Last Modified by GUIDE v2.5 06-May-2016 23:07:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SimApp_OpeningFcn, ...
                   'gui_OutputFcn',  @SimApp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SimApp is made visible.
function SimApp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimApp (see VARARGIN)

% Choose default command line output for SimApp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SimApp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SimApp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobuttonSine.
function radiobuttonSine_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonSine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonSine
value = get(hObject, 'value');
if value
    set(handles.radiobuttonGaussian, 'value', 0);
   
    set(handles.editFreqSine, 'Enable', 'on');
    set(handles.textFreqSine, 'Enable', 'on');
    set(handles.editONTimeSine, 'Enable', 'on');
    set(handles.textONTimeSine, 'Enable', 'on');
    
    set(handles.editFreqGaussian, 'Enable', 'off');
    set(handles.textFreqGaussian, 'Enable', 'off');
    set(handles.editBandwidthGaussian, 'Enable', 'off');
    set(handles.textBandwidthGaussian, 'Enable', 'off');
    
else
    set(handles.radiobuttonGaussian, 'value', 1);

    set(handles.editFreqSine, 'Enable', 'off');
    set(handles.textFreqSine, 'Enable', 'off');
    set(handles.editONTimeSine, 'Enable', 'off');
    set(handles.textONTimeSine, 'Enable', 'off');

    set(handles.editFreqGaussian, 'Enable', 'on');
    set(handles.textFreqGaussian, 'Enable', 'on');
    set(handles.editBandwidthGaussian, 'Enable', 'on');
    set(handles.textBandwidthGaussian, 'Enable', 'on');
end

% --- Executes on button press in radiobuttonGaussian.
function radiobuttonGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonGaussian
value = get(hObject, 'value');
if value
    set(handles.radiobuttonSine, 'value', 0);
  
    set(handles.editFreqGaussian, 'Enable', 'on');
    set(handles.textFreqGaussian, 'Enable', 'on');
    set(handles.editBandwidthGaussian, 'Enable', 'on');
    set(handles.textBandwidthGaussian, 'Enable', 'on');
    
    set(handles.editFreqSine, 'Enable', 'off');
    set(handles.textFreqSine, 'Enable', 'off');
    set(handles.editONTimeSine, 'Enable', 'off');
    set(handles.textONTimeSine, 'Enable', 'off');
    
else
    set(handles.radiobuttonSine, 'value', 1);
    
    set(handles.editFreqGaussian, 'Enable', 'off');
    set(handles.textFreqGaussian, 'Enable', 'off');
    set(handles.editBandwidthGaussian, 'Enable', 'off');
    set(handles.textBandwidthGaussian, 'Enable', 'off');
    
    set(handles.editFreqSine, 'Enable', 'on');
    set(handles.textFreqSine, 'Enable', 'on');
    set(handles.editONTimeSine, 'Enable', 'on');
    set(handles.textONTimeSine, 'Enable', 'on');
end


function editFreqSine_Callback(hObject, eventdata, handles)
% hObject    handle to editFreqSine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFreqSine as text
%        str2double(get(hObject,'String')) returns contents of editFreqSine as a double


% --- Executes during object creation, after setting all properties.
function editFreqSine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreqSine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFreqGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to editFreqGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFreqGaussian as text
%        str2double(get(hObject,'String')) returns contents of editFreqGaussian as a double


% --- Executes during object creation, after setting all properties.
function editFreqGaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreqGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBrowse1_Callback(hObject, eventdata, handles)
% hObject    handle to editBrowse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBrowse1 as text
%        str2double(get(hObject,'String')) returns contents of editBrowse1 as a double


% --- Executes during object creation, after setting all properties.
function editBrowse1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBrowse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'mytitle');
if filename ~= 0
    set(handles.editBrowse1, 'string', [pathname, filename]);
end

function editBrowse2_Callback(hObject, eventdata, handles)
% hObject    handle to editBrowse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBrowse2 as text
%        str2double(get(hObject,'String')) returns contents of editBrowse2 as a double


% --- Executes during object creation, after setting all properties.
function editBrowse2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBrowse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'mytitle');
if filename ~= 0
    set(handles.editBrowse2, 'string', [pathname, filename]);
end

function editONTimeSine_Callback(hObject, eventdata, handles)
% hObject    handle to editONTimeSine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editONTimeSine as text
%        str2double(get(hObject,'String')) returns contents of editONTimeSine as a double


% --- Executes during object creation, after setting all properties.
function editONTimeSine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editONTimeSine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBandwidthGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to editBandwidthGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBandwidthGaussian as text
%        str2double(get(hObject,'String')) returns contents of editBandwidthGaussian as a double


% --- Executes during object creation, after setting all properties.
function editBandwidthGaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBandwidthGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIterations_Callback(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIterations as text
%        str2double(get(hObject,'String')) returns contents of editIterations as a double


% --- Executes during object creation, after setting all properties.
function editIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGridSpacing_Callback(hObject, eventdata, handles)
% hObject    handle to editGridSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGridSpacing as text
%        str2double(get(hObject,'String')) returns contents of editGridSpacing as a double


% --- Executes during object creation, after setting all properties.
function editGridSpacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGridSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCourantFactor_Callback(hObject, eventdata, handles)
% hObject    handle to editCourantFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCourantFactor as text
%        str2double(get(hObject,'String')) returns contents of editCourantFactor as a double


% --- Executes during object creation, after setting all properties.
function editCourantFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCourantFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttonabc.
function radiobuttonabc_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonabc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonabc
value = get(hObject, 'value');
if value
    set(handles.editabcWidth, 'Enable', 'on');
    set(handles.textabcWidth, 'Enable', 'on');
    set(handles.editSigma0, 'Enable', 'on');
    set(handles.textSigma0, 'Enable', 'on');
else
    set(handles.editabcWidth, 'Enable', 'off');
    set(handles.textabcWidth, 'Enable', 'off');
    set(handles.editSigma0, 'Enable', 'off');
    set(handles.textSigma0, 'Enable', 'off');
end


function editabcWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editabcWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editabcWidth as text
%        str2double(get(hObject,'String')) returns contents of editabcWidth as a double


% --- Executes during object creation, after setting all properties.
function editabcWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editabcWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7



function editPlasmaFreq_Callback(hObject, eventdata, handles)
% hObject    handle to editPlasmaFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPlasmaFreq as text
%        str2double(get(hObject,'String')) returns contents of editPlasmaFreq as a double


% --- Executes during object creation, after setting all properties.
function editPlasmaFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPlasmaFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCollisionFreq_Callback(hObject, eventdata, handles)
% hObject    handle to editCollisionFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCollisionFreq as text
%        str2double(get(hObject,'String')) returns contents of editCollisionFreq as a double


% --- Executes during object creation, after setting all properties.
function editCollisionFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCollisionFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxSucceptibility_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxSucceptibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxSucceptibility as text
%        str2double(get(hObject,'String')) returns contents of editMaxSucceptibility as a double


% --- Executes during object creation, after setting all properties.
function editMaxSucceptibility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxSucceptibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxSourceField_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxSourceField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxSourceField as text
%        str2double(get(hObject,'String')) returns contents of editMaxSourceField as a double


% --- Executes during object creation, after setting all properties.
function editMaxSourceField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxSourceField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epsilon0 = (1e-9)/(36*pi);
mu0 = 4*pi*1e-7;

omegap = str2num(get(handles.editPlasmaFreq, 'string'));
vc = 2*pi*str2num(get(handles.editCollisionFreq, 'string'));
polarization = get(handles.popupmenu, 'value');

if polarization == 1
    TE = 1; TM = 0;
else
    TE = 0; TM = 1;
end

time_tot = str2num(get(handles.editIterations, 'string'));
N0 = str2num(get(handles.editConvolutionTerms, 'string'));
sigma0 = 1e5;
chimax = str2num(get(handles.editMaxSucceptibility, 'string'));

c = 3e8;
S = str2num(get(handles.editCourantFactor, 'string'));
delta = str2num(get(handles.editGridSpacing, 'string'));
deltat = (S/c)*delta;

sourceDrudePath = get(handles.editBrowse1, 'string');
I = imread(sourceDrudePath);
I = im2double(I);
[xdim, ydim] = size(I(:, :, 1)');
IB = I(:, :, 3);
IG = I(:, :, 2);
IR = I(:, :, 1);

drude = repmat((IB > 0.5), [1, 1, N0]);
sourceCoors = (IR > 0.5);

permittivityPath = get(handles.editBrowse2, 'string');
I1 = imread(permittivityPath);
I1 = im2double(I1);
IB1 = I1(:, :, 3);
IG1 = I1(:, :, 2);
IR1 = I1(:, :, 1);

epsx = epsilon0*(1 + IB1*chimax);
epsy = epsilon0*(1 + IR1*chimax);
epsz = epsilon0*(1 + IG1*chimax);

%convolution non-delta componnts
a = N0;
if N0 == 0
    a = 1;
end
epsxx = zeros(xdim, ydim, a);
epsyy = epsxx;
epszz = epsxx;

for i = 1:1:N0
  
    epsxx(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
    epsyy(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
    epszz(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
end

mux = mu0*ones(xdim, ydim);
muy = mu0*ones(xdim, ydim);
muz = mu0*ones(xdim, ydim);

sigma = zeros(xdim, ydim);

sourceWave = get(handles.radiobuttonSine, 'value');
if sourceWave == 1
    sine = 1; gaussian = 0;
    frequency = str2num(get(handles.editFreqSine, 'string'));
    onTime = str2num(get(handles.editONTimeSine, 'string'));
else
    sine = 0; gaussian = 1;
    centerFreq = str2num(get(handles.editFreqGaussian, 'string'));
    freqRatio = str2num(get(handles.editBandwidthGaussian, 'string'))/centerFreq;
end

zsource = zeros(1, time_tot);
maxSource = str2num(get(handles.editMaxSourceField, 'string'));
if gaussian == 1
    tc = gauspuls('cutoff', centerFreq, freqRatio,[],-60);
    t = -tc : deltat : tc;
    yi = gauspuls(t,centerFreq, freqRatio);
    sze = size(yi);
    
    for i = 1:1:sze(2)
           zsource(i) = maxSource*yi(i);
    end
    
elseif sine == 1
    for i = 1:1:onTime
        zsource(i) = maxSource*sin(2*pi*frequency*i*deltat);
    end
end

%ABC
ABCDecide = get(handles.radiobuttonabc, 'value');
if ABCDecide == 1
    
    x1 = str2num(get(handles.editabcWidth, 'string'));
    y1 = x1;
    x2 = xdim - x1; y2 = ydim - y1;

    ABC = ones(xdim, ydim);
    ABC(x1:x2, y1:y2) = 0;
    ABC_region = ABC > 0.5;
    
    sigma0 = str2num(get(handles.editSigma0, 'string'));

    for x = 1:1:x1
        sigma(x, y1:y2) = sigma0*(x/x1 - 1)^2;
    end
    for x = x2:1:xdim
        sigma(x, y1:y2) = sigma0*((x - x2)^2)/((x2-xdim)^2);
    end
    for y = 1:1:y1
        sigma(x1:x2, y) = sigma0*(y/y1 - 1)^2;
    end
    for y = y2:1:ydim
        sigma(x1:x2, y) = sigma0*((y - y2)^2)/((ydim-y2)^2);
    end

    sigma(1:x1, 1:y1) = (repmat(sigma(1:x1, y1+1), 1, y1) + repmat(sigma(x1+1, 1:y1), x1, 1));
    sigma(x2:xdim, 1:y1) = (repmat(sigma(x2:xdim, y1+1), 1, y1) + repmat(sigma(x1+1, 1:y1), xdim-x2+1, 1));
    sigma(x2:xdim, y2:ydim) = (repmat(sigma(x2:xdim, y1+1), 1, ydim-y2+1) + repmat(sigma(x1+1, y2:ydim), xdim-x2+1, 1));
    sigma(1:x1, y2:ydim) = (repmat(sigma(1:x1, y1+1), 1, ydim-y2+1) + repmat(sigma(x1+1, y2:ydim), x1, 1));

    abc = 1.05;

    epsx(ABC_region) = epsilon0*abc;
    epsy(ABC_region) = epsilon0*abc;
    epsz(ABC_region) = epsilon0*(1/abc);

    mux(ABC_region) = mu0*abc;
    muy(ABC_region) = mu0*abc;
    muz(ABC_region) = mu0*(1/abc);

    epsx(x2+1, y1:y2) = epsilon0;
    epsx(x1:x2, y2+1) = epsilon0;
    epsy(x2+1, y1:y2) = epsilon0;
    epsy(x1:x2, y2+1) = epsilon0;
    epsz(x2+1, y1:y2) = epsilon0;
    epsz(x1:x2, y2+1) = epsilon0;
end

Ez = zeros(xdim, ydim);
Ex = zeros(xdim, ydim);
Ey = zeros(xdim, ydim);
Hx = zeros(xdim, ydim);
Hy = zeros(xdim, ydim);
Hz = zeros(xdim, ydim);


figure;

plot(zsource);
title('source');
pause(1);

imagesc((sigma./sigma0)',[-1, 1]);
title('conductivity profile')
colorbar; %colormap(jet);\
pause(2);


imagesc((epszz(:, :, 1))',[-1, 1]);
title('drude material position')
colorbar; %colormap(jet);\
pause(2);


convStore = zeros(xdim, ydim, N0);
convStorex = zeros(xdim, ydim, N0);
convStorey = zeros(xdim, ydim, N0);

nx = 1; nxx = xdim-1;
ny = 1; nyy = ydim-1;

for n = 1:1:time_tot

    if TM == 1
        Hy(nx:nxx, ny:nyy) = Hy(nx:nxx,ny:nyy) + (0.5*deltat/delta)*(Ez(nx+1:nxx+1, ny:nyy) + Ez(nx+1:nxx+1, ny+1:nyy+1) - Ez(nx:nxx,ny:nyy) - Ez(nx:nxx, ny+1:nyy+1))./muy(nx:nxx,ny:nyy);
        Hx(nx:nxx, ny:nyy) = Hx(nx:nxx,ny:nyy) - (0.5*deltat/delta)*(Ez(nx:nxx, ny+1:nyy+1) + Ez(nx+1:nxx+1, ny+1:nyy+1) - Ez(nx:nxx,ny:nyy) - Ez(nx+1:nxx+1,ny:nyy))./mux(nx:nxx,ny:nyy);

        M = zeros(xdim-1, ydim-1);

        for i = 1:1:N0
            M = M + convStore(nx+1:nxx+1, ny+1:nyy+1, i).*epszz(nx+1:nxx+1, ny+1:nyy+1, i);    
        end

        Ez(nx+1:nxx+1, ny+1:nyy+1) = (1./(sigma(nx+1:nxx+1, ny+1:nyy+1)*deltat + epsz(nx+1:nxx+1, ny+1:nyy+1) + epszz(nx+1:nxx+1, ny+1:nyy+1, 1)*deltat*deltat)).*(epsz(nx+1:nxx+1, ny+1:nyy+1).*Ez(nx+1:nxx+1, ny+1:nyy+1) - deltat*deltat*M + (0.5*deltat/delta)*(Hy(nx+1:nxx+1, ny:nyy) + Hy(nx+1:nxx+1,ny+1:nyy+1) - Hy(nx:nxx,ny+1:nyy+1) - Hy(nx:nxx, ny:nyy) - Hx(nx+1:nxx+1,ny+1:nyy+1) - Hx(nx:nxx, ny+1:nyy+1) + Hx(nx+1:nxx+1,ny:nyy) + Hx(nx:nxx, ny:nyy)));

        convStore(:, :, 2:N0) = convStore(:, :, 1:N0-1);
        convStore(:, :, 1) = Ez;

        Ez(sourceCoors) = zsource(n);

        if mod(n, 2) == 1
            imagesc(Ez,[-1, 1]);
            colorbar; %colormap(jet);
        end
        getframe();
        
    else
       
        Mx = zeros(xdim-1, ydim-1);
        for i = 1:1:N0
            Mx = Mx + convStorex(nx+1:nxx+1, ny+1:nyy+1, i).*epsxx(nx+1:nxx+1, ny+1:nyy+1, i);
        end

        My = zeros(xdim-1, ydim-1);
        for i = 1:1:N0
            My = My + convStorey(nx+1:nxx+1, ny+1:nyy+1, i).*epsyy(nx+1:nxx+1, ny+1:nyy+1, i);
        end

        Ey(nx+1:nxx+1, ny+1:nyy+1) = (1./(epsy(nx+1:nxx+1, ny+1:nyy+1) + epsyy(nx+1:nxx+1, ny+1:nyy+1, 1)*deltat*deltat + sigma(nx+1: nxx+1, ny+1:nyy+1)*deltat)).*(epsy(nx+1:nxx+1,ny+1:nyy+1).*Ey(nx+1:nxx+1, ny+1:nyy+1) - My*deltat*deltat - (0.5*deltat/delta)*(Hz(nx+1:nxx+1, ny:nyy) + Hz(nx+1:nxx+1, ny+1:nyy+1) - Hz(nx:nxx,ny:nyy) - Hz(nx:nxx, ny+1:nyy+1)));
        Ex(nx+1:nxx+1, ny+1:nyy+1) = (1./(epsx(nx+1:nxx+1, ny+1:nyy+1) + epsxx(nx+1:nxx+1, ny+1:nyy+1, 1)*deltat*deltat + sigma(nx+1: nxx+1, ny+1:nyy+1)*deltat)).*(epsx(nx+1:nxx+1,ny+1:nyy+1).*Ex(nx+1:nxx+1, ny+1:nyy+1) - Mx*deltat*deltat + (0.5*deltat/delta)*(Hz(nx:nxx, ny+1:nyy+1) + Hz(nx+1:nxx+1, ny+1:nyy+1) - Hz(nx:nxx,ny:nyy) - Hz(nx+1:nxx+1, ny:nyy)));

        Hz(nx:nxx, ny:nyy) = Hz(nx:nxx, ny:nyy) + (0.5*deltat/delta)*(Ey(nx+1:nxx+1,ny+1:nyy+1) + Ey(nx+1:nxx+1, ny:nyy) - Ey(nx:nxx,ny+1:nyy+1)-Ey(nx:nxx,ny:nyy) - Ex(nx+1:nxx+1,ny+1:nyy+1) - Ex(nx:nxx, ny+1:nyy+1) + Ex(nx+1:nxx+1,ny:nyy) + Ex(nx:nxx, ny:nyy))./(-mux(nx:nxx, ny:nyy));

        convStorex(:, :, 2:N0) = convStorex(:, :, 1:N0-1);
        convStorex(:, :, 1) = Ex;

        convStorey(:, :, 2:N0) = convStorey(:, :, 1:N0-1);
        convStorey(:, :, 1) = Ey;

        Hz(sourceCoors) = zsource(n);

        if mod(n, 2) == 1
            imagesc(Hz,[-1, 1]);
            colorbar; %colormap(jet);
        end
        getframe();
        
    end
        
end

% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu


% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editConvolutionTerms_Callback(hObject, eventdata, handles)
% hObject    handle to editConvolutionTerms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editConvolutionTerms as text
%        str2double(get(hObject,'String')) returns contents of editConvolutionTerms as a double


% --- Executes during object creation, after setting all properties.
function editConvolutionTerms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConvolutionTerms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSigma0_Callback(hObject, eventdata, handles)
% hObject    handle to editSigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSigma0 as text
%        str2double(get(hObject,'String')) returns contents of editSigma0 as a double


% --- Executes during object creation, after setting all properties.
function editSigma0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
