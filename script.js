"use strict";
const N=6, names=["Tecnología","Operaciones","Marketing","Logística","Personal","Innovación"];
const truthIdeal=[1.2,0.8,1.5,1.1,0.9,1.3], truthStress=[1.0,0.7,1.25,1.3,0.8,1.1];
let A=[], b=[], x0=Array(N).fill(0), lastResult=null;
const $=id=>document.getElementById(id);
function matMul(A,x){return A.map(row=>row.reduce((s,v,j)=>s+v*x[j],0))}
function dot(a,b){return a.reduce((s,v,i)=>s+v*b[i],0)}
function norm(a){return Math.sqrt(dot(a,a))}
function sub(a,b){return a.map((v,i)=>v-b[i])}
function addScaled(a,b,t){return a.map((v,i)=>v+t*b[i])}
function scale(a,t){return a.map(v=>v*t)}
function transpose(A){return A[0].map((_,j)=>A.map(r=>r[j]))}
function norm2Matrix(A){return Math.max(...A.map(r=>r.reduce((s,v)=>s+Math.abs(v),0)))}
function identity(n){return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0))}
function makeC(diag=4,off=-1){let C=Array.from({length:N},(_,i)=>Array.from({length:N},(_,j)=>i===j?diag:(Math.abs(i-j)===1?off:0)));return C}
function makeScenario(kind){
 let C, d, truth;
 if(kind==="ideal"){C=makeC(4,-0.7);d=[1,1.1,.9,1.2,.95,1.05];truth=truthIdeal}
 else if(kind==="stress"){C=makeC(4.4,-1.25);d=[1,1.12,.94,1.18,1.02,1.1];truth=truthStress}
 else {C=makeC(4,-1);d=[1,3,10,30,100,300];truth=truthIdeal}
 const AA=C.map((row,i)=>row.map((v,j)=>v*d[i]*d[j]));
 return {A:AA,b:matMul(AA,truth),x0:Array(N).fill(0)};
}
function renderInputs(){
 $('matrixWrap').innerHTML='<table><thead><tr><th>A</th>'+names.map((_,i)=>`<th>x${i+1}</th>`).join('')+'</tr></thead><tbody>'+A.map((row,i)=>`<tr><td>${'x'+(i+1)}</td>${row.map((v,j)=>`<td><input type="number" step="any" data-a="${i},${j}" value="${fmtInput(v)}" aria-label="A ${i+1},${j+1}"></td>`).join('')}</tr>`).join('')+'</tbody></table>';
 $('bWrap').innerHTML='<table><thead><tr><th>Vector b</th><th>Valor</th></tr></thead><tbody>'+b.map((v,i)=>`<tr><td>b${i+1}</td><td><input type="number" step="any" data-b="${i}" value="${fmtInput(v)}"></td></tr>`).join('')+'</tbody></table>';
 $('x0Wrap').innerHTML='<table><thead><tr><th>Vector x₀</th><th>Valor inicial</th></tr></thead><tbody>'+x0.map((v,i)=>`<tr><td>x₀,${i+1}</td><td><input type="number" step="any" data-x0="${i}" value="${fmtInput(v)}"></td></tr>`).join('')+'</tbody></table>';
 document.querySelectorAll('[data-a]').forEach(el=>el.addEventListener('input',()=>{edited=true;let [i,j]=el.dataset.a.split(',').map(Number);A[i][j]=Number(el.value);if(i!==j&&$('keepSym').checked){A[j][i]=A[i][j];let other=document.querySelector(`[data-a="${j},${i}"]`);if(other)other.value=el.value;}runAll();}));
 document.querySelectorAll('[data-b]').forEach(el=>el.addEventListener('input',()=>{edited=true;b[Number(el.dataset.b)]=Number(el.value);runAll()}));
 document.querySelectorAll('[data-x0]').forEach(el=>el.addEventListener('input',()=>{edited=true;x0[Number(el.dataset.x0)]=Number(el.value);runAll()}));
}
function fmtInput(v){return Number(v.toPrecision(8)).toString()}
function fmt(v){if(!Number.isFinite(v))return '—';if(v!==0&&(Math.abs(v)>=1e5||Math.abs(v)<1e-4))return v.toExponential(4);return v.toFixed(6)}
function clearResults(){$('validation').textContent='Datos modificados. Pulsa “Validar y ejecutar” para recalcular.';$('validation').className='status';}
function loadScenario(kind){curKind=kind;edited=false;let s=makeScenario(kind);A=s.A.map(r=>r.slice());b=s.b.slice();x0=s.x0.slice();$('scenario').value=kind;renderInputs();clearResults();runAll()}
function isSymmetric(A,tol=1e-10){for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++)if(Math.abs(A[i][j]-A[j][i])>tol*Math.max(1,Math.abs(A[i][j]),Math.abs(A[j][i])))return false;return true}
function cholesky(A){let n=A.length,L=Array.from({length:n},()=>Array(n).fill(0));for(let i=0;i<n;i++)for(let j=0;j<=i;j++){let s=A[i][j];for(let k=0;k<j;k++)s-=L[i][k]*L[j][k];if(i===j){if(!(s>1e-12*Math.max(1,Math.abs(A[i][i]))))throw Error('La matriz no es definida positiva (falló Cholesky).');L[i][j]=Math.sqrt(s)}else L[i][j]=s/L[j][j]}return L}
function validateMatrix(A){if(!Array.isArray(A)||A.length===0||A.some(row=>!Array.isArray(row)||row.length!==A.length||row.some(v=>!Number.isFinite(v))))throw Error('A debe ser una matriz cuadrada con coeficientes numéricos finitos.');if(!isSymmetric(A))throw Error('La matriz A no es simétrica. Ajusta los términos aᵢⱼ y aⱼᵢ para que coincidan.');cholesky(A);return true}
function solveDirect(A,b){let n=A.length,M=A.map((r,i)=>r.concat([b[i]]));for(let k=0;k<n;k++){let p=k;for(let i=k+1;i<n;i++)if(Math.abs(M[i][k])>Math.abs(M[p][k]))p=i;if(Math.abs(M[p][k])<1e-15)throw Error('Matriz singular o casi singular.');[M[k],M[p]]=[M[p],M[k]];for(let i=k+1;i<n;i++){let f=M[i][k]/M[k][k];for(let j=k;j<=n;j++)M[i][j]-=f*M[k][j]}}let x=Array(n).fill(0);for(let i=n-1;i>=0;i--){let s=M[i][n];for(let j=i+1;j<n;j++)s-=M[i][j]*x[j];x[i]=s/M[i][i]}return x}
function eigenSym(A){let M=A.map(r=>r.slice()),n=M.length,maxIter=100*n*n;for(let it=0;it<maxIter;it++){let p=0,q=1,mx=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(M[i][j])>mx){mx=Math.abs(M[i][j]);p=i;q=j}if(mx<1e-14*Math.max(1,...M.map((r,i)=>Math.abs(r[i]))))break;let tau=(M[q][q]-M[p][p])/(2*M[p][q]);let t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau));let c=1/Math.sqrt(1+t*t),s=t*c,app=M[p][p],aqq=M[q][q],apq=M[p][q];M[p][p]=app-t*apq;M[q][q]=aqq+t*apq;M[p][q]=M[q][p]=0;for(let k=0;k<n;k++)if(k!==p&&k!==q){let mkp=M[k][p],mkq=M[k][q];M[k][p]=M[p][k]=c*mkp-s*mkq;M[k][q]=M[q][k]=s*mkp+c*mkq}}return M.map((r,i)=>r[i])}
function condSPD(A){let ev=eigenSym(A),lo=Math.min(...ev),hi=Math.max(...ev);return lo>0?hi/lo:Infinity}
function transformedJacobi(A){let d=A.map((r,i)=>r[i]);return A.map((row,i)=>row.map((v,j)=>v/Math.sqrt(d[i]*d[j])))}
function cg(A,b,xinit,tol,maxit,precond){let x=xinit.slice(),r=sub(b,matMul(A,x)),bn=norm(b)||1,history=[norm(r)/bn],iterations=0;let z=precond?r.map((v,i)=>v/A[i][i]):r.slice(),p=z.slice(),rz=dot(r,z),converged=history[0]<tol;
 for(let k=0;k<maxit&&!converged;k++){let Ap=matMul(A,p),den=dot(p,Ap);if(!(den>0)||!Number.isFinite(den))break;let alpha=rz/den;x=addScaled(x,p,alpha);r=addScaled(r,Ap,-alpha);let rel=norm(r)/bn;history.push(rel);iterations=k+1;if(rel<tol){converged=true;break}z=precond?r.map((v,i)=>v/A[i][i]):r.slice();let rzNew=dot(r,z);if(!Number.isFinite(rzNew)||Math.abs(rz)<1e-300)break;let beta=rzNew/rz;p=addScaled(z,p,beta);rz=rzNew}
 let trueR=norm(sub(b,matMul(A,x)))/bn;return {x,iterations,history,relativeResidual:trueR,converged:trueR<tol}}
function makeTable(headers,rows){return '<table><thead><tr>'+headers.map(h=>`<th>${h}</th>`).join('')+'</tr></thead><tbody>'+rows.map(row=>'<tr>'+row.map(v=>`<td>${v}</td>`).join('')+'</tr>').join('')+'</tbody></table>'}
function runAll(){
 try{
 if(!b.every(Number.isFinite)||!x0.every(Number.isFinite))throw Error('Los vectores b y x₀ deben contener números finitos.');
 validateMatrix(A);
 let tol=Number($('tol').value),maxit=Number($('maxit').value);if(!(tol>0&&tol<1))throw Error('La tolerancia debe estar entre 0 y 1.');if(!(maxit>=1&&Number.isInteger(maxit)))throw Error('El máximo de iteraciones debe ser un entero positivo.');
 let ref=solveDirect(A,b),simple=cg(A,b,x0,tol,maxit,false),usePre=$('preconditioner').value==='jacobi',pcg=cg(A,b,x0,tol,maxit,usePre);
 let kA=condSPD(A),kP=usePre?condSPD(transformedJacobi(A)):kA;
 lastResult={ref,simple,pcg,kA,kP,usePre};
 $('validation').className='status ok';$('validation').textContent='Validación correcta: A es simétrica y Cholesky confirma que es definida positiva.';
 $('kA').textContent=fmt(kA);$('kP').textContent=fmt(kP);$('iters').textContent=`${simple.iterations} / ${pcg.iterations}`;$('resid').textContent=fmt(pcg.relativeResidual);
 $('resultsStatus').className='status '+(pcg.converged?'ok':'warn');$('resultsStatus').textContent=`GC: ${simple.converged?'convergió':'no alcanzó tolerancia'}; GCP: ${pcg.converged?'convergió':'no alcanzó tolerancia'}. Método directo calculado por eliminación de Gauss con pivoteo parcial.`;
 $('solutionTable').innerHTML=makeTable(['Variable','Referencia directa','GC','GCP'],names.map((n,i)=>[`${n} (x${i+1})`,fmt(ref[i]),fmt(simple.x[i]),fmt(pcg.x[i])]));
 $('errorTable').innerHTML=makeTable(['Variable','Referencia','GCP','Ea','Er','Ep (%)'],names.map((n,i)=>{let ea=Math.abs(ref[i]-pcg.x[i]),er=ref[i]===0?null:ea/Math.abs(ref[i]);return [`${n} (x${i+1})`,fmt(ref[i]),fmt(pcg.x[i]),fmt(ea),er===null?'No definido':fmt(er),er===null?'No definido':fmt(er*100)]}));
 $('comparison').innerHTML=makeTable(['Indicador','Ideal','Bajo estrés','Mal condicionado'],buildComparison());
 drawChart(simple.history,pcg.history);
 renderExtras(ref,pcg,usePre);renderConclusions(simple,pcg,kA,kP);
 }catch(e){clearStale();$('validation').className='status err';$('validation').textContent=e.message;$('resultsStatus').className='status err';$('resultsStatus').textContent='No se calcularon resultados: corrige los datos y vuelve a ejecutar.'}
}
function buildComparison(){let scenarios=['ideal','stress','ill'],out=[];let keys=['Solución referencia','Solución GCP','Error absoluto GCP','Error porcentual GCP','Iteraciones GC','Iteraciones GCP','κ₂(A)','κ₂ precondicionada','Residuo relativo GCP','¿Convergió GCP?','Interpretación'];let results=scenarios.map(kind=>{let s=makeScenario(kind),ref=solveDirect(s.A,s.b),gc=cg(s.A,s.b,s.x0,1e-6,100,false),pc=cg(s.A,s.b,s.x0,1e-6,100,true);return {ref,gc,pc,k:condSPD(s.A),kp:condSPD(transformedJacobi(s.A)),kind}});
for(let key of keys){let vals=results.map(r=>{if(key==='Solución referencia')return '['+r.ref.map(fmt).join(', ')+']';if(key==='Solución GCP')return '['+r.pc.x.map(fmt).join(', ')+']';if(key==='Error absoluto GCP')return '['+r.ref.map((v,i)=>fmt(Math.abs(v-r.pc.x[i]))).join(', ')+']';if(key==='Error porcentual GCP')return '['+r.ref.map((v,i)=>v===0?'ND':fmt(Math.abs(v-r.pc.x[i])/Math.abs(v)*100)).join(', ')+']';if(key==='Iteraciones GC')return r.gc.iterations;if(key==='Iteraciones GCP')return r.pc.iterations;if(key==='κ₂(A)')return fmt(r.k);if(key==='κ₂ precondicionada')return fmt(r.kp);if(key==='Residuo relativo GCP')return fmt(r.pc.relativeResidual);if(key==='¿Convergió GCP?')return r.pc.converged?'Sí':'No';return TXT[r.kind]});out.push([key,...vals])}return out}
function drawChart(h1,h2){let c=$('chart'),ctx=c.getContext('2d'),w=c.width,h=c.height;ctx.clearRect(0,0,w,h);ctx.fillStyle='#0c1626';ctx.fillRect(0,0,w,h);let pad={l:65,r:20,t:22,b:48},pw=w-pad.l-pad.r,ph=h-pad.t-pad.b;let vals=h1.concat(h2).filter(v=>v>0&&Number.isFinite(v));let maxY=Math.max(0,...vals.map(Math.log10)),minY=Math.min(-8,...vals.map(Math.log10));if(maxY-minY<1)minY=maxY-1;ctx.strokeStyle='#304158';ctx.fillStyle='#a7b6cc';ctx.font='13px system-ui';for(let i=0;i<=5;i++){let y=pad.t+ph*i/5;ctx.beginPath();ctx.moveTo(pad.l,y);ctx.lineTo(w-pad.r,y);ctx.stroke();let v=maxY-(maxY-minY)*i/5;ctx.fillText('10^'+v.toFixed(1),5,y+4)}let maxX=Math.max(1,h1.length-1,h2.length-1);for(let i=0;i<=maxX;i++){let x=pad.l+pw*i/maxX;if(i%Math.max(1,Math.ceil(maxX/8))===0){ctx.fillText(String(i),x-4,h-20)}}ctx.fillText('Iteración',w/2-25,h-5);function line(arr,color){ctx.strokeStyle=color;ctx.lineWidth=2.5;ctx.beginPath();arr.forEach((v,i)=>{let lv=Math.log10(Math.max(v,1e-16)),x=pad.l+pw*i/maxX,y=pad.t+ph*(maxY-lv)/(maxY-minY);if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)});ctx.stroke()}line(h1,'#8bb7ff');line(h2,'#66d9c2');ctx.fillStyle='#8bb7ff';ctx.fillRect(w-195,12,12,4);ctx.fillStyle='#e8eef8';ctx.fillText('GC',w-177,19);ctx.fillStyle='#66d9c2';ctx.fillRect(w-125,12,12,4);ctx.fillStyle='#e8eef8';ctx.fillText('GCP',w-107,19)}

let curKind='ideal',edited=false;
const TXT={
ideal:"Escenario estable: costos propios moderados y sinergias débiles entre áreas. El equilibrio asigna 1,2; 0,8; 1,5; 1,1; 0,9 y 1,3 MM Bs.; un cambio pequeño en una meta apenas se propaga a las demás y el cálculo es rápido y exacto (κ bajo).",
stress:"La inflación (+10 % en costos propios) y la mayor interdependencia, junto con la caída de ingresos, obligan a recortar la asignación de la mayoría de las áreas; Logística sube por costos de transporte. El condicionamiento casi no cambia: el estrés desplaza el punto de equilibrio, no vuelve difícil el problema numérico.",
ill:"Las áreas se miden en escalas monetarias muy distintas: κ(A) ≈ 10⁵ aunque la economía del modelo sea la del caso ideal. Es un problema de escalamiento, no de riesgo financiero; Jacobi lo corrige (κ baja a ≈ 2,6) y recupera las mismas asignaciones del escenario ideal."};
function minorsOf(A){const L=cholesky(A);let p=1;return L.map((r,i)=>{p*=r[i]*r[i];return p})}
function scenInfo(kind){const s=makeScenario(kind),ev=eigenSym(s.A);return{s,lo:Math.min(...ev),hi:Math.max(...ev),mn:minorsOf(s.A),ref:solveDirect(s.A,s.b),k:condSPD(s.A)}}
function clearStale(){['kA','kP','iters','resid','totalBudget','valorOrg'].forEach(id=>$(id).textContent='—');['solutionTable','errorTable','refTable','cashTable','rulesTable','spdCurrent','sysEq'].forEach(id=>$(id).innerHTML='');$('spdCurrent').className='status err';$('conclusionText').innerHTML='<p>No hay resultados: corrige los datos para generar conclusiones.</p>';let c=$('chart');c.getContext('2d').clearRect(0,0,c.width,c.height)}
function renderExtras(ref,pcg,usePre){
 $('sysEq').textContent=A.map((r,i)=>r.map((v,j)=>({v,j})).filter(t=>t.v!==0).map((t,k)=>(k?(t.v<0?' − ':' + '):(t.v<0?'−':''))+(+Math.abs(t.v).toPrecision(5))+'·x'+(t.j+1)).join('')+'  =  '+(+b[i].toPrecision(5))).join('\n');
 const kinds=['ideal','stress','ill'],lab={ideal:'Ideal',stress:'Bajo estrés',ill:'Mal condicionado'},info=kinds.map(scenInfo);
 $('spdTable').innerHTML=makeTable(['Escenario','A = Aᵀ','λmín','λmáx','Menores principales (Sylvester)','¿Definida positiva?'],info.map((o,i)=>[lab[kinds[i]],isSymmetric(o.s.A)?'Sí':'No',fmt(o.lo),fmt(o.hi),o.mn.map(fmt).join(' · '),o.lo>0&&o.mn.every(v=>v>0)?'Sí ✔':'No ✘']));
 const ev=eigenSym(A),mn=minorsOf(A);$('spdCurrent').className='status ok';$('spdCurrent').textContent=`Datos cargados ahora: A simétrica; λmín = ${fmt(Math.min(...ev))} > 0; menores principales todos > 0 (${mn.map(fmt).join(' · ')}). Por Sylvester, A es definida positiva.`;
 $('refTable').innerHTML=makeTable(['Variable','Área','x de referencia (MM Bs.)'],names.map((n,i)=>[`x${i+1}`,n,fmt(ref[i])]))+`<p class="small muted">Residuo relativo de la referencia: ${(norm(sub(b,matMul(A,ref)))/(norm(b)||1)).toExponential(2)}. Cholesky se usa solo para verificar la definición positiva.</p>`;
 $('mWrap').textContent=usePre?`Jacobi: M = diag(A) = diag(${A.map((r,i)=>fmt(r[i])).join(', ')}). En cada paso Mz = r se resuelve dividiendo cada componente de r por a_ii.`:'Sin precondicionador: M = I (GC simple).';
 const tot=ref.reduce((s,v)=>s+v,0)||1,neg=ref.filter(v=>v<0).length;
 $('totalBudget').textContent=fmt(tot)+' MM Bs.';$('valorOrg').textContent=fmt(-0.5*dot(b,ref));
 $('cashTable').innerHTML=makeTable(['Área','Asignación x (MM Bs.)','% del total'],names.map((n,i)=>[n,fmt(ref[i]),fmt(ref[i]/tot*100)+' %']));
 $('rulesTable').innerHTML=makeTable(['Regla','Descripción','Estado'],[['R1 · Asignaciones no negativas','x_i ≥ 0 en todas las áreas',neg?`✘ ${neg} área(s) con asignación negativa: replantea A o b`:'✔ Cumple'],['R2 · Simetría de interdependencias','a_ij = a_ji (la dependencia entre áreas es mutua)','✔ Cumple (verificado)'],['R3 · Estabilidad del costo','A definida positiva: el equilibrio es un mínimo único','✔ Cumple (Cholesky + Sylvester)'],['R4 · Precisión','Residuo relativo < tolerancia',pcg.converged?'✔ Cumple':'✘ No cumple']]);
 const a=info[0],c=info[1],ratio=c.k/a.k,ta=a.ref.reduce((s,v)=>s+v,0),tc=c.ref.reduce((s,v)=>s+v,0);
 $('stressTable').innerHTML=makeTable(['Área','x Ideal','x Bajo estrés','Variación'],names.map((n,i)=>[n,fmt(a.ref[i]),fmt(c.ref[i]),((c.ref[i]/a.ref[i]-1)*100).toFixed(1)+' %']).concat([['Total Σx',fmt(ta),fmt(tc),((tc/ta-1)*100).toFixed(1)+' %']]));
 $('stressText').innerHTML=`κ(A) pasa de <b>${fmt(a.k)}</b> a <b>${fmt(c.k)}</b> (×${ratio.toFixed(2)}). ${ratio<2?'El estrés no altera significativamente el condicionamiento: ambos valores son bajos y el algoritmo converge igual de rápido.':'El estrés sí deteriora el condicionamiento.'} El efecto principal es financiero: el equilibrio se desplaza y el presupuesto total cambia ${((tc/ta-1)*100).toFixed(1)} %.`}
function renderConclusions(simple,pcg,kA,kP){
 $('conclusionText').innerHTML=generalConclusions()+`<h3>Resultado de los datos cargados actualmente</h3><ul class="list"><li><b>Técnica.</b> κ(A) = ${fmt(kA)} y κ(M⁻¹A) = ${fmt(kP)}. GC usó ${simple.iterations} iteraciones y GCP ${pcg.iterations}; residuo relativo GCP = ${fmt(pcg.relativeResidual)} (${pcg.converged?'convergió':'no convergió'}). Con n = 6, GC termina en a lo sumo 6 pasos en aritmética exacta, por lo que la ventaja de Jacobi se ve sobre todo en el número de condición y en los casos mal escalados.</li><li><b>Financiera.</b> ${edited?'Datos editados: compara las asignaciones con las metas y revisa las reglas de negocio de la sección 3.':TXT[curKind]}</li></ul><p class="small muted">Los resultados son demostrativos (datos sintéticos) y no constituyen una recomendación de inversión.</p>`}

function generalConclusions(){
 const S=['ideal','stress','ill'].map(k=>{const s=makeScenario(k),ref=solveDirect(s.A,s.b),g=cg(s.A,s.b,s.x0,1e-6,100,false),p=cg(s.A,s.b,s.x0,1e-6,100,true);return{ref,g,p,kA:condSPD(s.A),kP:condSPD(transformedJacobi(s.A)),tot:ref.reduce((x,v)=>x+v,0),val:-0.5*dot(s.b,ref),ea:Math.max(...ref.map((v,i)=>Math.abs(v-p.x[i])))}});
 const [a,c,m]=S,e=v=>v.toExponential(1),pc=v=>((v)*100).toFixed(1).replace('.',',')+' %',k=v=>v>=1000?String(Math.round(v)):v.toFixed(2).replace('.',','),
 maxEa=Math.max(...S.map(o=>o.ea)),maxRes=Math.max(...S.map(o=>o.p.relativeResidual)),ch=c.ref.map((v,i)=>v/a.ref[i]-1),low=ch.map((v,i)=>v<=Math.min(...ch)+1e-6?i:-1).filter(i=>i>=0),f2=v=>v.toFixed(2).replace('.',','),iMax=ch.indexOf(Math.max(...ch));
 return `<h3>Conclusiones generales de los tres escenarios</h3><ul class="list">
<li><b>Exactitud.</b> En los tres escenarios GCP converge y reproduce la solución de referencia (método directo): error absoluto máximo ${e(maxEa)} y residuo relativo máximo ${e(maxRes)}, muy por debajo de la tolerancia 10⁻⁶. El método es preciso y confiable para sistemas simétricos definidos positivos como este.</li>
<li><b>Condicionamiento.</b> κ(A) = ${k(a.kA)} (ideal), ${k(c.kA)} (bajo estrés) y ${k(m.kA)} (mal condicionado). El estrés multiplica κ por solo ${(c.kA/a.kA).toFixed(2).replace('.',',')}, así que no vuelve difícil el problema: únicamente desplaza el punto de equilibrio. El mal escalamiento, en cambio, lo multiplica por ≈ ${Math.round(m.kA/a.kA).toLocaleString('es')}, y ese deterioro es numérico, no financiero: viene de mezclar unidades de magnitudes distintas.</li>
<li><b>Velocidad de convergencia.</b> Iteraciones GC / GCP: ${a.g.iterations}/${a.p.iterations} (ideal), ${c.g.iterations}/${c.p.iterations} (estrés) y ${m.g.iterations}/${m.p.iterations} (mal condicionado). Jacobi reduce κ de ${k(m.kA)} a ${k(m.kP)}, del mismo orden que en el caso ideal. Con n = 6 la ganancia en iteraciones es pequeña porque GC termina en a lo sumo 6 pasos; en sistemas más grandes las iteraciones crecen aproximadamente con √κ, y ahí la reducción de κ sí se traduce en una ventaja grande.</li>
<li><b>Interpretación financiera del estrés.</b> Frente al caso ideal, la inflación de costos, la mayor interdependencia y la caída de ingresos reducen el presupuesto total del equilibrio de ${f2(a.tot)} a ${f2(c.tot)} MM Bs. (${pc(c.tot/a.tot-1)}) y el valor organizacional f(x*) pasa de ${f2(a.val)} a ${f2(c.val)} (menos negativo = menos valor neto; ${pc(c.val/a.val-1)}). El modelo recomienda recortar más en ${low.map(i=>names[i]).join(' y ')} (${pc(Math.min(...ch))} cada una) y reforzar ${names[iMax]} (${pc(ch[iMax])}), cuyo costo de transporte sube. Es decir: la empresa debe reajustar su plan, pero el cálculo sigue siendo estable.</li>
<li><b>Interpretación del mal condicionamiento.</b> Las asignaciones de equilibrio son las mismas que en el caso ideal; lo único que cambia es la escala de cada variable. La lección para un sistema real es normalizar las unidades antes de modelar o usar un precondicionador diagonal, que lo hace automáticamente.</li>
<li><b>Aporte de Ingeniería de Sistemas.</b> El modelo traduce procesos interdepartamentales en coeficientes de A, y las reglas de negocio (asignaciones no negativas, simetría, estabilidad del costo y precisión) se validan automáticamente en cada ejecución.</li>
<li><b>Alcance y limitaciones.</b> Los datos son sintéticos, el modelo es cuadrático y n = 6; los resultados ilustran el método y no constituyen una recomendación de inversión.</li></ul>`}
function solveAll(AA,bb,xx,tol){const ref=solveDirect(AA,bb),g=cg(AA,bb,xx,tol,100,false),p=cg(AA,bb,xx,tol,100,true);return{ref,g,p,Ea:Math.max(...ref.map((v,i)=>Math.abs(v-p.x[i])))}}
function runQA(){
 const rows=[],tol=1e-6,add=(c,e,o,ok)=>rows.push([c,e,o,ok?'<b class="good">Correcto</b>':'<b class="bad">Fallido</b>']);
 for(const [kind,lab] of [['ideal','Ideal'],['stress','Bajo estrés'],['ill','Mal condicionado']]){const s=makeScenario(kind);let msg,ok=false;try{validateMatrix(s.A);const r=solveAll(s.A,s.b,s.x0,tol);ok=r.p.converged&&r.Ea<1e-6;msg=`A simétrica definida positiva; GCP convergió en ${r.p.iterations} iteraciones; Ea máx. = ${r.Ea.toExponential(2)}`}catch(e){msg=e.message}add(`1. Cargar escenario ${lab}`,'A válida, GCP converge y coincide con la referencia (Ea < 10⁻⁶)',msg,ok)}
 const t=(M,p)=>{try{validateMatrix(M);return['No se detectó el error',false]}catch(e){return[e.message,e.message.includes(p)]}};
 let ns=identity(N);ns[0][1]=2;let r=t(ns,'no es simétrica');add('2a. Matriz no simétrica (a₁₂ = 2, a₂₁ = 0)','Mensaje de error sin bloquear la página',r[0],r[1]);
 let np=identity(N);np[0][1]=np[1][0]=2;r=t(np,'definida positiva');add('2b. Simétrica no definida positiva (a₁₂ = a₂₁ = 2, diagonal 1)','Mensaje de error sin bloquear la página',r[0],r[1]);
 let nn=identity(N);nn[2][2]=NaN;r=t(nn,'finitos');add('2c. Coeficiente no numérico (NaN)','Mensaje de error sin bloquear la página',r[0],r[1]);
 {const s=makeScenario('ideal'),base=solveAll(s.A,s.b,s.x0,tol),A2=s.A.map(q=>q.slice()),b2=s.b.slice();A2[0][0]+=1;b2[0]+=0.1;let msg,ok=false;try{validateMatrix(A2);const q=solveAll(A2,b2,s.x0,tol);ok=Math.abs(q.ref[0]-base.ref[0])>1e-6&&q.p.converged&&q.Ea<1e-6&&q.p.relativeResidual<tol;msg=`x₁ de referencia: ${fmt(base.ref[0])} → ${fmt(q.ref[0])}; GCP la reproduce (Ea máx. = ${q.Ea.toExponential(2)}, residuo = ${q.p.relativeResidual.toExponential(2)})`}catch(e){msg=e.message}add('3. Edición dinámica: a₁₁ + 1 y b₁ + 0,1 en el escenario Ideal','Se recalculan referencia, GCP y errores; la solución cambia y GCP la sigue',msg,ok)}
 {const s=makeScenario('ideal'),q=solveAll(s.A,s.b,[5,-3,2,0,1,4],tol);add('3b. Vector inicial x₀ = (5, −3, 2, 0, 1, 4)','GCP converge desde un x₀ distinto de cero',`GCP: ${q.p.iterations} iteraciones; Ea máx. = ${q.Ea.toExponential(2)}`,q.p.converged&&q.Ea<1e-6)}
 {const s=makeScenario('ill'),q=solveAll(s.A,s.b,s.x0,tol);add('4. Precondicionamiento en el escenario Mal condicionado','Iteraciones GCP ≤ iteraciones GC',`GC = ${q.g.iterations}; GCP = ${q.p.iterations}; κ(A) = ${fmt(condSPD(s.A))}; κ(M⁻¹A) = ${fmt(condSPD(transformedJacobi(s.A)))}`,q.p.iterations<=q.g.iterations&&q.p.converged)}
 $('qaTable').innerHTML=makeTable(['Entrada / caso de prueba','Resultado esperado','Resultado obtenido','Conclusión'],rows)}
function goResults(){const v=$('validation');if(v.className.includes('ok')){const t=$('resultados');t.scrollIntoView({behavior:'smooth',block:'start'});t.classList.remove('flash');void t.offsetWidth;t.classList.add('flash')}else v.scrollIntoView({behavior:'smooth',block:'center'})}
document.querySelectorAll('[data-load-scenario]').forEach(btn=>btn.addEventListener('click',()=>{loadScenario(btn.dataset.loadScenario);goResults()}));
$('loadScenario').addEventListener('click',()=>{loadScenario($('scenario').value);goResults()});$('run').addEventListener('click',()=>{runAll();goResults()});$('reset').addEventListener('click',()=>loadScenario($('scenario').value));$('scenario').addEventListener('change',()=>loadScenario($('scenario').value));['tol','maxit','preconditioner'].forEach(id=>$(id).addEventListener('change',runAll));
$('runQA').addEventListener('click',runQA);
loadScenario('ideal');runQA();


