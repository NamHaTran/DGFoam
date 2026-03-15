Bản backup 22/05/2021:
- Đổi hoàn toàn sang mô hình Extended NSF Durst.
- Cần modify lại điều kiện biên vận tốc tại wall: vận tốc convection vẫn như classical NSF, vận tốc diffusion pháp tuyến với wall bằng 0.
- MÔ HÌNH DURST CHẠY OK KHÔNG BỊ CRASHHHHHHHH! REESE-GREENSHIELD LÀ HÀNG ĐỂUUUUUU.

Bản backup 27/05/2021:
- Update mô hình Durst:
	+ Dễ cài đặt vào/gỡ ra khỏi Classical NSF-DG, code của Classical NSF chỉ bị can thiệp vào ở 3 hàm.
	+ Bổ sung tích phân bị thiếu ở RHS khi giải ở bậc cao (ngoài additional term ở viscous component còn tích phân int(div(Ec)*phi) trên volume Omega. Vì div(Ec) ~ sum(Ec_i*div(phi_i)) nên khi order i=1 thì term này bằng 0) => hơi tốn não khi tìm cách tính tích phân này.
	+ Điều chỉnh điều kiện biên của vận tốc tại wall (bằng hàm bcForExtNSF_Durst::dropNormDiffVel cũng thuộc mô hình Durst) để loại bỏ thành phần vận tốc do self diffusion vuông góc với wall (đương nhiên là hiện tượng self diffusion không xảy ra trên phương vuông góc wall).
	
Bản backup 28/05/2021:
- Cải thiện hiệu suất của mô hình Durst: nhân hệ số theta2 của positivity preserving limiter vào command tính Ec_order trong hàm correctEnergyEqnVolIntTerm() vì ở bậc cao, Ec_order có thể có giá trị unphysical ở vùng có strong discontinuity -> làm solver mất ổn định.
- Sau khi cải thiện, solver có thể chạy mô hình Durst Kn=0.01 ở số Co = 0.1 (lúc trước chỉ chạy được ở Co < 0.03).

Bản backup 31/05/2021:
- Sửa bug lưu mDx và mDy vào array bị sai vị trí. Bug này làm case Kn=0.05 và 0.1 bị sai nặng vì ở vùng số Kn này ảnh hưởng của diffusive flux mD khá lớn.

Bản backup 14/06/2021:
- Add thêm điều kiện biên zeroGradRhoUncorrectP cho file p, điều kiện biên này apply rho- = rho+ và grad(rho)- = grad(rho)+, không correct p theo phương trình khí lý tưởng. Điều kiện biên này giúp tính grad(rho) tại biên wall tốt hơn (khí bị nén ở sát wall nên gradient theo chiều hướng về wall, còn với điều kiện biên interpFrmDensity, chiều gradient sát wall bị ngược lại).
- Add thêm hệ số Dm = 1/Sc vào thành phần self diffusion (1.32 đối với Argon). Điều này làm ảnh hưởng của self diffusion mạnh hơn. Nếu dùng như original formula của tài liệu cho mô hình Durst, effect của self diffusion không đủ mạnh để tạo ra khác biệt giữa mô hình Classical và Durst.

Bản backup 24/08/2021:
Mô hình Durst
- Chốt điều kiện biên tại wall trong trường hợp self-diffusion
	+ Vận tốc diffusion (gồm vận tốc do mass diffusion và thermophoresis (Sorret term)) theo phương pháp tuyến bằng 0.
	+ Vận tốc diffusion do thermophoresis bằng 0 theo cả pháp tuyến và tiếp tuyến bền mặt (nên trong code, dTx và dTy set bằng 0 tại biên).
	+ Điều kiện của T và rho giống classical NSF.

Điều kiện biên
- Add folder điều kiện biên custom, các điều kiện biên mới sau này có thể được thêm vào dưới dạng điều kiện biên custom.
- Add điều kiện biên interior:
	+ Với các field T, U: giá trị scalar và gradient phía - bằng giá trị phía +
	+ Với field rho: giá trị scalar phía - bằng phía +. Giá trị gradient phía - thu được khi apply hàm zero normal gradient cho grad phía +.
- Add bộ điều kiện biên non-equilibrium, hiện tại gồm điều kiện Maxwell slip và Smoluchowsky T Jump. Các điều kiện biên non-equilibrium sau này có thể được add vào tương tự như 2 đk biên trước.

Fix bug
- Fix bug giải kết quả 'nan' của hàm tìm hình chiếu vuông góc của cell centroid lên 1 cạnh.

- Code update lên github cùng ngày.

Bản backup 14/09/2021:
Điều kiện biên
- Update điều kiện biên Maxwell và Smoluchowsky:
	+ Đổi cách lưu surfaceField data, mảng lưu giá trị T và U tại surface giờ là nonEqmSurfaceField, lưu giá trị tại từng điểm Gauss trên edge.
	+ Giá trị T và U tính từ đk biên nonequilibrium tại từng edge giờ là giá trị tại từng điểm Gauss, chứ không phải là giá trị tại hình chiếu vuông góc của centroid xuống edge nữa.
	+ Update hàm giải đk biên. Trong đó đạo hàm tại các điểm Gauss được tính bằng (phi_Gauss - phi_C)*(dot(i,n))/delta, với
		* phi là biến đang xét. phi_Gauss là biến đang xét tại điểm Gauss trên surface, có thể là biến hoặc giá trị có sẵn.
		* dot(i,n) là dot product của i (vector đơn vị chỉ phương của đoạn thẳng C_GaussPt) và n (vector pháp tuyến đơn vị của edge).
		* delta là độ dài đoạn C_GaussPt.
	Phương pháp này giúp bậc chính xác của kết quả tính tăng lên theo bậc của bài toán và số điểm Gauss, thay vì chỉ là bậc 1 như trước (giá trị T, U tại các điểm Gauss được coi như bằng giá trị T, U tại hình chiếu của centroid xuống cạnh).
	+ Update hàm đọc file *surface.txt, có thể bỏ qua nếu file không tồn tại, hoặc đọc kết quả từ file theo format cũ rồi tính giá trị trung bình và phân phối lại giá trị này tới các điểm Gauss trên cạnh. Update này giúp có thể chạy tiếp từ kết quả tính từ version cũ của code mà không bị conflict.
	+ Update thêm 1 số hàm xử lý hình học mới cần cho nonequilibrium BC.
Fix bug
- Fix bug nhỏ về chính tả trong message output ra terminal.

- Code update lên github cùng ngày.


Bản backup 29/09/2021:
Điều kiện biên
- Update điều kiện biên Smoluchowsky:
	+ Thêm hàm SmoluchowskyTJump::calcTJump_FDMTypeSemiImplicit để giải điều kiện biên Smoluchowsky theo cách semi-implicit. Trong đó hệ số nhớt được tính từ TJump cũ bằng hàm tính hệ số nhớt, vì vậy hàm semi-implicit này có thể dùng cho mọi mô hình nhớt.
	+ Hàm SmoluchowskyTJump::calcTJump_FDMTypeSemiImplicit là hàm được dùng trong hàm nonEquilibriumBCs::updateBCs thay cho hàm calcTJump_FDMTypeImplicit (gây sai số khi dùng với mô hình nhớt khác Sutherland).
- Update điều kiện biên của mô hình ENSE DURST:
	+ Tại biên wall, sau khi tính các term grad(T) và grad(rho), phải gán grad(T)=0 (BỎ SORRET TERM trong thành phần self-diffusion tại wall). 
	+ Sau khi tính self-diffusion flux tại phía Minus của wall, phải reflect vector flux này lại bằng hàm bcForExtNSF_Durst::removeNormTermOfVector để remove thành phần self-diffusion flux vuông góc wall.
	+ Chú ý rằng phải reflect self-diffusion flux chứ không phải reflect grad(rho) hay grad(T). Luôn reflect self-diffusion flux trong mọi trường hợp, kể cả trường hợp self-diffusion flux theo phương vuông góc wall hướng ra khỏi wall (dot product của self-diffusion flux là vector pháp tuyến n < 0).
	+ Cách apply điều kiện biên tại wall này đã được test trên case cylinder Kn=0.05 và cho kết quả tốt, đồng thời case khá ổn định và có thể chạy với Max Co = 0.2.
- Add điều kiện biên zeroRhoGrad để triệt tiêu normal grad(rho) tại wall (không quan trọng).

- Code update lên github cùng ngày.


Bản backup 23/11/2021:
- Thêm scheme Euler cho time discretization. Với case cylinder Kn=0.1, chạy bằng Euler nhanh hơn, ổn định hơn và có thể dùng được HLLE flux.
- Fix bug nhỏ về chính tả trong message output ra terminal.

- Code update lên github cùng ngày.


Bản backup 12/7/2022:
- Tách numberOfGaussPoints thành numberOfGaussPoints1D và numberOfGaussPoints2D.
- Thêm nhóm điều kiện biên timeVaryingBCs để chứa các điều kiện biên thay đổi theo thời gian. Chuyển field lưu giá trị điều kiện biên của điều kiện nonEquilibriumBCs từ nonEqmSurfaceField sang SurfaceBCFields để thống nhất với các điều kiện biên thuộc nhóm timeVaryingBCs trong tương lai.
- Thêm điều kiện biên waveTransmissive thuộc timeVaryingBCs. Hiện tại chạy OK với vector field U, nhưng bị crash khi chạy với scalar field.