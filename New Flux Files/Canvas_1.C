void Canvas_1()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Wed Jan 22 00:21:48 2020) by ROOT version 6.18/04
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",522,189,2034,1125);
   Canvas_1->Range(-15,-5.568162e-11,135,5.011345e-10);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Double_t xAxis1[961] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.125, 4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5, 5.125, 5.25, 5.375, 5.5, 5.625, 5.75, 5.875, 6, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8, 8.125, 8.25, 8.375, 8.5, 8.625, 8.75, 8.875, 9, 9.125, 9.25, 9.375, 9.5, 9.625, 9.75, 9.875, 10, 10.125, 10.25, 10.375, 10.5, 10.625, 10.75, 10.875, 11, 11.125, 11.25, 11.375, 11.5, 11.625, 11.75, 11.875, 12, 12.125, 12.25, 12.375, 12.5, 12.625, 12.75, 12.875, 13, 13.125, 13.25, 13.375, 13.5, 13.625, 13.75, 13.875, 14, 14.125, 14.25, 14.375, 14.5, 14.625, 14.75, 14.875, 15, 15.125, 15.25, 15.375, 15.5, 15.625, 15.75, 15.875, 16, 16.125, 16.25, 16.375, 16.5, 16.625, 16.75, 16.875, 17, 17.125, 17.25, 17.375, 17.5, 17.625, 17.75, 17.875, 18, 18.125, 18.25, 18.375, 18.5, 18.625, 18.75, 18.875, 19, 19.125, 19.25, 19.375, 19.5, 19.625, 19.75, 19.875, 20, 20.125, 20.25, 20.375, 20.5, 20.625, 20.75, 20.875, 21, 21.125, 21.25, 21.375, 21.5, 21.625, 21.75, 21.875, 22, 22.125, 22.25, 22.375, 22.5, 22.625, 22.75, 22.875, 23, 23.125, 23.25, 23.375, 23.5, 23.625, 23.75, 23.875, 24, 24.125, 24.25, 24.375, 24.5, 24.625, 24.75, 24.875, 25, 25.125, 25.25, 25.375, 25.5, 25.625, 25.75, 25.875, 26, 26.125, 26.25, 26.375, 26.5, 26.625, 26.75, 26.875, 27, 27.125, 27.25, 27.375, 27.5, 27.625, 27.75, 27.875, 28, 28.125, 28.25, 28.375, 28.5, 28.625, 28.75, 28.875, 29, 29.125, 29.25, 29.375, 29.5, 29.625, 29.75, 29.875, 30, 30.125, 30.25, 30.375, 30.5, 30.625, 30.75, 30.875, 31, 31.125, 31.25, 31.375, 31.5, 31.625, 31.75, 31.875, 32, 32.125, 32.25, 32.375, 32.5, 32.625, 32.75, 32.875, 33, 33.125, 33.25, 33.375, 33.5, 33.625, 33.75, 33.875, 34, 34.125, 34.25, 34.375, 34.5, 34.625, 34.75, 34.875, 35, 35.125, 35.25, 35.375, 35.5, 35.625, 35.75, 35.875, 36, 36.125, 36.25, 36.375, 36.5, 36.625, 36.75, 36.875, 37, 37.125, 37.25, 37.375, 37.5, 37.625, 37.75, 37.875, 38, 38.125, 38.25, 38.375, 38.5, 38.625, 38.75, 38.875, 39, 39.125, 39.25, 39.375, 39.5, 39.625, 39.75, 39.875, 40, 40.125, 40.25, 40.375, 40.5, 40.625, 40.75, 40.875, 41, 41.125, 41.25, 41.375, 41.5, 41.625, 41.75, 41.875, 42, 42.125, 42.25, 42.375, 42.5, 42.625, 42.75, 42.875, 43, 43.125, 43.25, 43.375, 43.5, 43.625, 43.75, 43.875, 44, 44.125, 44.25, 44.375, 44.5, 44.625, 44.75, 44.875, 45, 45.125, 45.25, 45.375, 45.5, 45.625, 45.75, 45.875, 46, 46.125, 46.25, 46.375, 46.5, 46.625, 46.75, 46.875, 47, 47.125, 47.25, 47.375, 47.5, 47.625, 47.75, 47.875, 48, 48.125, 48.25, 48.375, 48.5, 48.625, 48.75, 48.875, 49, 49.125, 49.25, 49.375, 49.5, 49.625, 49.75, 49.875, 50, 50.125, 50.25, 50.375, 50.5, 50.625, 50.75, 50.875, 51, 51.125, 51.25, 51.375, 51.5, 51.625, 51.75, 51.875, 52, 52.125, 52.25, 52.375, 52.5, 52.625, 52.75, 52.875, 53, 53.125, 53.25, 53.375, 53.5, 53.625, 53.75, 53.875, 54, 54.125, 54.25, 54.375, 54.5, 54.625, 54.75, 54.875, 55, 55.125, 55.25, 55.375, 55.5, 55.625, 55.75, 55.875, 56, 56.125, 56.25, 56.375, 56.5, 56.625, 56.75, 56.875, 57, 57.125, 57.25, 57.375, 57.5, 57.625, 57.75, 57.875, 58, 58.125, 58.25, 58.375, 58.5, 58.625, 58.75, 58.875, 59, 59.125, 59.25, 59.375, 59.5, 59.625, 59.75, 59.875, 60, 60.125, 60.25, 60.375, 60.5, 60.625, 60.75, 60.875, 61, 61.125, 61.25, 61.375, 61.5, 61.625, 61.75, 61.875, 62, 62.125, 62.25, 62.375, 62.5, 62.625, 62.75, 62.875, 63, 63.125, 63.25, 63.375, 63.5, 63.625, 63.75, 63.875, 64, 64.125, 64.25, 64.375, 64.5, 64.625, 64.75, 64.875, 65, 65.125, 65.25, 65.375, 65.5, 65.625, 65.75, 65.875, 66, 66.125, 66.25, 66.375, 66.5, 66.625, 66.75, 66.875, 67, 67.125, 67.25, 67.375, 67.5, 67.625, 67.75, 67.875, 68, 68.125, 68.25, 68.375, 68.5, 68.625, 68.75, 68.875, 69, 69.125, 69.25, 69.375, 69.5, 69.625, 69.75, 69.875, 70, 70.125, 70.25, 70.375, 70.5, 70.625, 70.75, 70.875, 71, 71.125, 71.25, 71.375, 71.5, 71.625, 71.75, 71.875, 72, 72.125, 72.25, 72.375, 72.5, 72.625, 72.75, 72.875, 73, 73.125, 73.25, 73.375, 73.5, 73.625, 73.75, 73.875, 74, 74.125, 74.25, 74.375, 74.5, 74.625, 74.75, 74.875, 75, 75.125, 75.25, 75.375, 75.5, 75.625, 75.75, 75.875, 76, 76.125, 76.25, 76.375, 76.5, 76.625, 76.75, 76.875, 77, 77.125, 77.25, 77.375, 77.5, 77.625, 77.75, 77.875, 78, 78.125, 78.25, 78.375, 78.5, 78.625, 78.75, 78.875, 79, 79.125, 79.25, 79.375, 79.5, 79.625, 79.75, 79.875, 80, 80.125, 80.25, 80.375, 80.5, 80.625, 80.75, 80.875, 81, 81.125, 81.25, 81.375, 81.5, 81.625, 81.75, 81.875, 82, 82.125, 82.25, 82.375, 82.5, 82.625, 82.75, 82.875, 83, 83.125, 83.25, 83.375, 83.5, 83.625, 83.75, 83.875, 84, 84.125, 84.25, 84.375, 84.5, 84.625, 84.75, 84.875, 85, 85.125, 85.25, 85.375, 85.5, 85.625, 85.75, 85.875, 86, 86.125, 86.25, 86.375, 86.5, 86.625, 86.75, 86.875, 87, 87.125, 87.25, 87.375, 87.5, 87.625, 87.75, 87.875, 88, 88.125, 88.25, 88.375, 88.5, 88.625, 88.75, 88.875, 89, 89.125, 89.25, 89.375, 89.5, 89.625, 89.75, 89.875, 90, 90.125, 90.25, 90.375, 90.5, 90.625, 90.75, 90.875, 91, 91.125, 91.25, 91.375, 91.5, 91.625, 91.75, 91.875, 92, 92.125, 92.25, 92.375, 92.5, 92.625, 92.75, 92.875, 93, 93.125, 93.25, 93.375, 93.5, 93.625, 93.75, 93.875, 94, 94.125, 94.25, 94.375, 94.5, 94.625, 94.75, 94.875, 95, 95.125, 95.25, 95.375, 95.5, 95.625, 95.75, 95.875, 96, 96.125, 96.25, 96.375, 96.5, 96.625, 96.75, 96.875, 97, 97.125, 97.25, 97.375, 97.5, 97.625, 97.75, 97.875, 98, 98.125, 98.25, 98.375, 98.5, 98.625, 98.75, 98.875, 99, 99.125, 99.25, 99.375, 99.5, 99.625, 99.75, 99.875, 100, 100.125, 100.25, 100.375, 100.5, 100.625, 100.75, 100.875, 101, 101.125, 101.25, 101.375, 101.5, 101.625, 101.75, 101.875, 102, 102.125, 102.25, 102.375, 102.5, 102.625, 102.75, 102.875, 103, 103.125, 103.25, 103.375, 103.5, 103.625, 103.75, 103.875, 104, 104.125, 104.25, 104.375, 104.5, 104.625, 104.75, 104.875, 105, 105.125, 105.25, 105.375, 105.5, 105.625, 105.75, 105.875, 106, 106.125, 106.25, 106.375, 106.5, 106.625, 106.75, 106.875, 107, 107.125, 107.25, 107.375, 107.5, 107.625, 107.75, 107.875, 108, 108.125, 108.25, 108.375, 108.5, 108.625, 108.75, 108.875, 109, 109.125, 109.25, 109.375, 109.5, 109.625, 109.75, 109.875, 110, 110.125, 110.25, 110.375, 110.5, 110.625, 110.75, 110.875, 111, 111.125, 111.25, 111.375, 111.5, 111.625, 111.75, 111.875, 112, 112.125, 112.25, 112.375, 112.5, 112.625, 112.75, 112.875, 113, 113.125, 113.25, 113.375, 113.5, 113.625, 113.75, 113.875, 114, 114.125, 114.25, 114.375, 114.5, 114.625, 114.75, 114.875, 115, 115.125, 115.25, 115.375, 115.5, 115.625, 115.75, 115.875, 116, 116.125, 116.25, 116.375, 116.5, 116.625, 116.75, 116.875, 117, 117.125, 117.25, 117.375, 117.5, 117.625, 117.75, 117.875, 118, 118.125, 118.25, 118.375, 118.5, 118.625, 118.75, 118.875, 119, 119.125, 119.25, 119.375, 119.5, 119.625, 119.75, 119.875, 120}; 
   
   TH1D *enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1 = new TH1D("enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1","",960, xAxis1);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(1,2.919944e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(2,6.69766e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(3,7.545952e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(4,9.761691e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(5,1.037763e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(6,1.245707e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(7,1.494912e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(8,1.68101e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(9,1.989947e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(10,1.948708e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(11,2.171184e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(12,2.383776e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(13,2.683703e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(14,3.050558e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(15,3.202379e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(16,3.074057e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(17,3.192083e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(18,3.314705e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(19,3.823715e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(20,3.616416e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(21,3.907589e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(22,3.967539e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(23,3.943239e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(24,4.043252e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(25,4.050619e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(26,3.949709e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(27,4.242409e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(28,4.004748e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(29,3.777892e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(30,3.790272e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(31,4.037349e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(32,3.915102e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(33,3.668317e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(34,3.67365e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(35,3.69982e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(36,3.395406e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(37,3.358796e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(38,3.356348e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(39,3.290598e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(40,2.964463e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(41,2.965483e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(42,3.019951e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(43,2.822597e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(44,2.734641e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(45,2.431299e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(46,2.554558e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(47,2.392048e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(48,2.366451e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(49,2.069351e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(50,1.914422e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(51,2.004668e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(52,1.880752e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(53,1.731294e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(54,1.780565e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(55,1.788873e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(56,1.380845e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(57,1.422805e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(58,1.370298e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(59,1.425452e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(60,1.194184e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(61,1.364021e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(62,1.161713e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(63,1.179552e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(64,1.036957e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(65,1.095215e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(66,1.048989e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(67,1.057942e-10);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(68,8.254731e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(69,9.337477e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(70,7.537972e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(71,9.065161e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(72,8.32103e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(73,7.687828e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(74,8.668751e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(75,8.058622e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(76,7.817044e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(77,6.940856e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(78,6.863172e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(79,6.693765e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(80,7.542561e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(81,6.721434e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(82,6.15027e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(83,7.118429e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(84,5.902703e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(85,6.395405e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(86,6.261992e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(87,6.377255e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(88,5.683439e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(89,5.618385e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(90,6.128012e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(91,5.258721e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(92,5.361488e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(93,4.893686e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(94,4.772387e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(95,5.224846e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(96,4.535117e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(97,5.222207e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(98,5.053267e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(99,3.920373e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(100,4.862875e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(101,5.725119e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(102,4.766287e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(103,4.503846e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(104,4.71025e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(105,4.381742e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(106,3.968659e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(107,4.016386e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(108,4.214698e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(109,4.037267e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(110,3.999514e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(111,3.722257e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(112,4.091281e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(113,3.206748e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(114,3.752016e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(115,4.046759e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(116,3.095586e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(117,3.290374e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(118,3.413352e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(119,3.115395e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(120,3.490157e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(121,3.230761e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(122,2.884766e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(123,3.069597e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(124,2.829425e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(125,2.893683e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(126,2.610111e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(127,3.578462e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(128,2.920962e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(129,2.171518e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(130,2.119908e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(131,2.832824e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(132,2.544814e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(133,3.498482e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(134,2.129386e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(135,1.905638e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(136,2.026165e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(137,2.129547e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(138,2.475874e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(139,2.80176e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(140,2.220479e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(141,1.776391e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(142,1.757115e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(143,2.453085e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(144,2.068437e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(145,2.135088e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(146,1.679856e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(147,1.842114e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(148,1.829275e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(149,1.481626e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(150,1.711932e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(151,1.309187e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(152,1.631591e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(153,1.399597e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(154,1.454198e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(155,1.332664e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(156,1.330657e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(157,1.325124e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(158,1.606824e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(159,1.225422e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(160,1.419078e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(161,1.467132e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(162,1.082755e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(163,1.273179e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(164,1.476645e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(165,1.207626e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(166,1.816559e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(167,9.626574e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(168,1.308997e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(169,1.544094e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(170,1.064834e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(171,1.075347e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(172,9.565088e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(173,8.267366e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(174,1.130568e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(175,1.12195e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(176,1.077218e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(177,1.438203e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(178,9.970317e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(179,8.884846e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(180,7.17149e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(181,8.545807e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(182,8.183427e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(183,1.168703e-11);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(184,7.019809e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(185,8.024519e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(186,8.94447e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(187,6.142077e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(188,4.821276e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(189,7.806726e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(190,8.077926e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(191,7.798104e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(192,6.498064e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(193,8.44806e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(194,4.60042e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(195,7.795415e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(196,7.254641e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(197,6.240015e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(198,9.004507e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(199,5.659043e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(200,4.41643e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(201,5.809431e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(202,8.467773e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(203,5.31602e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(204,6.439218e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(205,3.123265e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(206,7.044707e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(207,5.272803e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(208,3.375995e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(209,3.860476e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(210,4.993396e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(211,4.498151e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(212,5.021704e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(213,3.556446e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(214,7.861011e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(215,4.10395e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(216,4.736925e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(217,5.393282e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(218,2.008713e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(219,4.756533e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(220,2.250922e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(221,2.411316e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(222,1.62707e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(223,5.619681e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(224,3.228321e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(225,1.844356e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(226,5.429518e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(227,3.477183e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(228,2.354473e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(229,2.795449e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(230,3.028891e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(231,3.850352e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(232,5.355834e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(233,6.35057e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(234,4.574875e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(235,4.604154e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(236,3.433288e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(237,3.405232e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(238,1.509695e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(239,3.684263e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(240,3.475293e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(241,3.516392e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(242,1.70678e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(243,2.515517e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(244,4.944888e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(245,1.61534e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(246,3.222359e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(247,2.562633e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(248,1.801916e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(249,1.806727e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(250,3.175057e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(251,2.326686e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(252,1.894472e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(253,3.51593e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(254,4.589566e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(255,2.126586e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(256,1.568006e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(257,2.732937e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(258,1.936409e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(259,2.365514e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(260,2.210321e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(261,4.304753e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(262,2.960363e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(263,1.741766e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(264,1.744408e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(265,4.421377e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(266,1.109806e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(267,2.279484e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(268,1.710451e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(269,1.257639e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(270,2.071013e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(271,1.808334e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(272,1.593151e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(273,1.190996e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(274,4.467965e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(275,2.625648e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(276,2.485512e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(277,1.446776e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(278,2.766667e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(279,2.138601e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(280,4.840673e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(282,1.799655e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(283,1.634906e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(284,2.615847e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(285,1.679421e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(287,4.76797e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(288,1.261295e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(290,1.277335e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(291,6.636846e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(292,1.714959e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(293,1.979669e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(295,5.546938e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(296,1.04517e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(297,8.804788e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(298,6.445383e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(299,6.00782e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(302,4.716224e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(303,1.118311e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(304,5.461087e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(305,1.982373e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(306,1.42305e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(307,1.053927e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(308,1.414721e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(309,1.318589e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(311,9.39508e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(314,8.450306e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(316,5.548544e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(318,8.502617e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(319,6.145862e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(320,1.133274e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(322,5.220115e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(326,1.182225e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(327,4.333874e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(328,7.78062e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(331,8.505362e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(332,7.029423e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(333,1.00557e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(334,5.668104e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(337,6.325669e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(342,2.307699e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(344,6.302129e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(347,8.846371e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(349,8.839274e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(350,1.385041e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(351,9.156346e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(354,9.792044e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(356,2.114031e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(371,1.185454e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(373,2.103808e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(383,1.012365e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(384,1.053893e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(390,9.174707e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(395,9.857565e-13);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(399,1.395553e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(412,1.520738e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(421,1.441884e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(441,1.693963e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(456,1.033084e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(487,1.813313e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetBinContent(535,1.652942e-12);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetEntries(2179907);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 2179907");
   ptstats_LaTex = ptstats->AddText("Mean  =  6.235");
   ptstats_LaTex = ptstats->AddText("Std Dev   =  5.586");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->SetLineColor(ci);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetTitle("Energy(GeV) 960 bins for 0-120GeV");
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetLabelFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetLabelSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetTitleSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetTitleOffset(1);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetXaxis()->SetTitleFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetYaxis()->SetLabelFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetYaxis()->SetLabelSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetYaxis()->SetTitleSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetYaxis()->SetTitleOffset(1);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetYaxis()->SetTitleFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetZaxis()->SetLabelFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetZaxis()->SetLabelSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetZaxis()->SetTitleSize(0.035);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetZaxis()->SetTitleOffset(1);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->GetZaxis()->SetTitleFont(42);
   enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS__1->Draw("");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
