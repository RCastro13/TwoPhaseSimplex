#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <fstream>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

const double EPS = 1e-9; //constante para evitar erro numérico

void printTable(const Matrix& table){
    for(const auto& row : table){
        for(const auto& val : row){
            cout << setw(10) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

//função para verificar se a solução é ótima
bool isOptimal(const Matrix& table){
    //olho se todos os elementos do meu custo são < -EPS
    const Vector& objective = table.back();
    for(int i=0;i<objective.size()-1;i++){
        if(objective[i] < -EPS){
            return false;
        }
    }
    return true;
}

//função para selecionar a coluna pivô
int getPivotColumn(const Matrix& table){
    const Vector& objective = table.back();
    int pivotCol = 0;
    //seleciono o menor elemento do meu custo
    for(int i = 1; i < objective.size() - 1; ++i){
        if(objective[i] < objective[pivotCol]){
            pivotCol = i;
        }
    }
    return pivotCol;
}

//função para selecionar a linha pivô
int getPivotRow(const Matrix& table, int pivotCol){
    int pivotRow = -1;
    double minRatio = numeric_limits<double>::infinity();
    //caso eu não encontre linha válida retorno -1
    for(int i = 0; i < table.size() - 1; ++i){
        double element = table[i][pivotCol];
        //considera-se somente razões com elementos positivos
        if(element > EPS){ 
            double ratio = table[i].back() / element;
            if(ratio >= 0 && ratio < minRatio){
                minRatio = ratio;
                pivotRow = i;
            }
        }
    }
    return pivotRow;
}

//função para pivotear a linha
void pivot(Matrix& table, int pivotRow, int pivotCol, Matrix& optimalAuxiliar){
    double pivotElement = table[pivotRow][pivotCol];

    //divisão da linha pivô pelo elemento pivô
    for(int i=0;i<table[0].size();i++){
        if(table[pivotRow][i] != 0){
            table[pivotRow][i] = table[pivotRow][i] / pivotElement;
        }
    }
    //replicação da operação na matriz do certificado
    for(int i=0;i<optimalAuxiliar[0].size();i++){
        if(optimalAuxiliar[pivotRow][i] != 0){
            optimalAuxiliar[pivotRow][i] = optimalAuxiliar[pivotRow][i] / pivotElement;
        }
    }

    //atualizando as outras linhas para zerar a coluna do pivô
    for(int i = 0; i < table.size(); ++i){
        if(i != pivotRow){
            double factor = table[i][pivotCol];
            for(int j = 0; j < table[i].size(); ++j){
                // double calculus = factor * table[pivotRow][j];
                // if(calculus != 0) 
                table[i][j] -= factor * table[pivotRow][j];
            }
            //aplicando mesma operação na matriz do certificado
            for(int h=0;h<optimalAuxiliar[0].size();++h){
                // double calculus = factor * optimalAuxiliar[pivotRow][h];
                // if(calculus != 0) 
                optimalAuxiliar[i][h] -= factor * optimalAuxiliar[pivotRow][h];
            }
        }
    }
}

//função para checar se minha PL é ilimitada
bool checkUnbounded(Matrix& table, Vector& unboudedCertificate, int pivotCol){
    bool isUnbounded = true;
    //caso a minha coluna pivô seja toda negativa, a PL é ilimitada
    for(int i=0;i<table.size()-1;i++){
        if(table[i][pivotCol] > EPS){
            isUnbounded = false;
            return isUnbounded;
        }
        else{
            continue;
        }
    }

    return isUnbounded;
}

//função para construir meu certificado de ilimitado
void buildUnbounded(Matrix& table, Vector& unboudedCertificate, int pivotCol){
    unboudedCertificate[pivotCol] = 1;
    int rows = table.size() - 1;

    //dado que estou no meu formato canônico, eu consigo montar o certificado identificando colunas básicas e não básicas
    for(int j = 0; j < table[0].size() - 1; ++j){
        int zeros = 0;
        int uns = 0;
        
        for(int i=0;i<table.size()-1;i++){
            if(table[i][j] == 1) uns++;
            else if(table[i][j] == 0) zeros++;
        }
        if(uns ==1 && zeros==table.size()-2 && j!=pivotCol){
            
            int pos = 0;
            for(int h=0;h<table.size()-1;h++){
                if(table[h][j] == 1){
                    pos = h;
                    break;
                }
            }
            unboudedCertificate[j] = -1*table[pos][pivotCol];
        }
    }

}

//função principal do Simplex
int Simplex(Matrix& table, Matrix& optimalAuxiliar, int viability, int flag, vector<int>& zeroColumns, Vector& unboudedCertificate){
    int counter = 0;
    vector<int> baseLines;
    //resolvendo a auxiliar - 1ª FASE
    if(flag==0){
        while(!isOptimal(table)){
            int pivotCol = getPivotColumn(table);
            int pivotRow = getPivotRow(table, pivotCol);

            pivot(table, pivotRow, pivotCol, optimalAuxiliar);
        }

        if(table[table.size()-1][table[0].size()-1] < -EPS){
            // cout << "Problema original inviavel" << endl;
        }
        else{
            // cout << "Solucao otima encontrada." << endl;
            viability = 1;
        }
    }
    //resolvendo a original - 2ª FASE
    else{
        while(!isOptimal(table)){
            viability = 1;
            if(counter == 0){
                counter = 1;
                Vector pivotsRow(zeroColumns.size());
                for(int i=0;i<zeroColumns.size();i++){
                    pivotsRow[i] = i;
                }
                //canonizando a original a partir das bases da auxiliar
                
                while(zeroColumns.size()!=0){
                    int pivotCol = zeroColumns[0];
                    int pivotRow = getPivotRow(table, pivotCol);
                    
                    // int pivotRow = pivotsRow[0];
                    if(pivotRow == -1){
                        zeroColumns.erase(zeroColumns.begin());
                        pivotsRow.erase(pivotsRow.begin());
                        continue;
                    }
                    
                    zeroColumns.erase(zeroColumns.begin());
                    pivotsRow.erase(pivotsRow.begin());
                    
                    if(checkUnbounded(table, unboudedCertificate, pivotCol)){
                        buildUnbounded(table, unboudedCertificate, pivotCol);
                        viability = 0;
                        return viability;
                    }
                    
                    pivot(table, pivotRow, pivotCol, optimalAuxiliar);
                    
                }
            }
            else{
                int pivotCol = getPivotColumn(table);
                int pivotRow = getPivotRow(table, pivotCol);
            
                if(pivotRow == -1){
                    if(checkUnbounded(table, unboudedCertificate, pivotCol)){
                        buildUnbounded(table, unboudedCertificate, pivotCol);
                        viability = 0;
                        return viability;
                    }
                    // continue;
                }
                pivot(table, pivotRow, pivotCol, optimalAuxiliar);
            }
        }
    }
    return viability;
}

//função para canonizar minha matriz auxiliar
void canonize(Matrix& table){
    int numRows = table.size();
    int numCols = table[0].size();

    //ajuste da linha do custo subtraindo as linhas de restrição correspondentes
    for(int i = 0; i < numRows - 1; ++i){
        for(int j = 0; j < numCols; ++j){
            table[numRows - 1][j] -= table[i][j];
        }
    }
}

//função para identificar as colunas básicas encontradas pela auxiliar
vector<int> findZeroColumns(const Matrix& table, int n, int m){
    int numRows = table.size();
    int numCols = table[0].size();
    vector<int> zeroColumns;
    vector<bool> cols(n+m, false);
    
    //analisar a última linha, exceto o valor da função
    for(int j = 0; j < numCols - 1; ++j){
        int zeros = 0;
        int uns = 0;
        if(table[numRows - 1][j] == 0){
            for(int i=0;i<table.size()-1;i++){
                if(table[i][j] == 1) uns++;
                else if(table[i][j] == 0) zeros++;
            }
            if(uns==1 && zeros==table.size()-2){
                zeroColumns.push_back(j);
                cols[j] = true;
            }
        }
    }

    //caso existe uma coluna referente a uma variável artifical eu troco por uma variável original
    for(int i=0;i<zeroColumns.size();i++){
        if(zeroColumns[i] > m-1){
            for(int j=0;j<cols.size()-n;j++){
                if(cols[j] == false){
                    zeroColumns[i] = j;
                    cols[j] = true;
                    cols[i] = false;
                }
            }
        }
    }

    return zeroColumns;
}

//função para identificar os valores de cada variável na solução
Vector findIdentityColumnsValues(const Matrix& table){
    int numRows = table.size();
    int numCols = table[0].size();
    vector<double> auxiliaryVector(numCols - 1, 0.0);

    //identifico onde estão as colunas da identidade e pega o valor de custo e atribuo em um vetor na posição respectiva
    for(int j = 0; j < numCols - 1; ++j){
        bool isIdentityColumn = true;
        int identityRow = -1;

        for(int i = 0; i < numRows - 1; ++i){
            if(table[i][j] == 1){
                if(identityRow == -1){
                    identityRow = i;
                } else{
                    isIdentityColumn = false;
                    break;
                }
            } else if (table[i][j] != 0){
                isIdentityColumn = false;
                break;
            }
        }

        //se for uma coluna da identidade eu atualizo o vetor auxiliar
        if(isIdentityColumn && identityRow != -1){
            auxiliaryVector[j] = table[identityRow][numCols - 1];
        }
    }

    return auxiliaryVector;
}

int main(int argc, char* argv[]){
    //entrada por argumento
    if(argc < 2){
        cerr << "Uso: " << argv[0] << " <arquivo_entrada>" << endl;
        return 1;
    }
    string filename = argv[1];
    ifstream inputFile(filename);
    if(!inputFile){
        cerr << "Erro ao abrir o arquivo " << filename << endl;
        return 1;
    }

    //entrada
    int n, m;
    inputFile >> n >> m;
    
    Vector cost(m);
    for(int i = 0; i < m; ++i){
        inputFile >> cost[i];
    }

    Matrix tableVarAuxiliares(n + 1, Vector(m + n + 1, 0.0));
    
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            inputFile >> tableVarAuxiliares[i][j];
        }
        inputFile >> tableVarAuxiliares[i][m + n];
        //variáveis de folga:
        tableVarAuxiliares[i][m + i] = 1.0;
    }

    inputFile.close();

    //*-1 o vetor de custo para o tableau
    for(int j = 0; j < m; ++j){
        tableVarAuxiliares[n][j] = -cost[j];
    }

    //montando a tabela do simplex da original
    Matrix table(n + 1, Vector(m + 1, 0.0));
    for(int i=0; i<n+1;i++){
        for(int j=0;j<m;j++){
            table[i][j] = tableVarAuxiliares[i][j]; 
        }
        table[i][m] = tableVarAuxiliares[i][m+n];
    }

    //matriz que vai guardar os certificados de ótimo
    Matrix optimalCertificate(n + 1, Vector(n, 0.0));
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(i == j){
                optimalCertificate[i][j] = 1;
            }
        }
    }

    Matrix optimalAuxiliar = optimalCertificate;

    //multiplicando a linha por -1 caso algum b negativo
    for(int i = 0; i < n; ++i){
        if(tableVarAuxiliares[i][m + n] < 0){
            for(int j = 0; j < m; ++j){
                if(tableVarAuxiliares[i][j] != 0){
                    tableVarAuxiliares[i][j] = tableVarAuxiliares[i][j] * -1;
                    //aplicando mesma operação na matriz do certificado
                    for(int h=0;h<n;++h){
                        if(optimalAuxiliar[i][h] != 0){
                            optimalAuxiliar[i][h] = optimalAuxiliar[i][h] * -1;
                        }
                    }
                }
            }
            tableVarAuxiliares[i][n+m] = tableVarAuxiliares[i][n+m] * -1;
        }
    }

    //montando a tabela do simplex da auxiliar
    Matrix auxiliarTable = tableVarAuxiliares;
    for(int i=0;i<m+n+1;i++){
        if(i < m || i==m+n){
            auxiliarTable[n][i] = 0;
        }else if(i>=m && i!=m+n){
            auxiliarTable[n][i] = 1;
        }
    }
    canonize(auxiliarTable);
    canonize(optimalAuxiliar);

    int viability = 0; //variável de controle da minha resposta
    int isAuxiliarOrOriginal = 0; //0=auxiliar, 1=original
    Vector unboudedCertificate(m, 0.0); //certificado de ilimitado

    //resolução da PL auxiliar
    vector<int> zeroColumns = findZeroColumns(auxiliarTable, n, m);
    int result = Simplex(auxiliarTable, optimalAuxiliar, viability, isAuxiliarOrOriginal, zeroColumns, unboudedCertificate);
    
    //result = 1 (PL auxiliar com valor ótimo 0), = 0 (PL auxiliar com valor ótimo < 0))
    if(result){
        isAuxiliarOrOriginal = 1;
        zeroColumns = findZeroColumns(auxiliarTable, n, m);
        //resolução da PL original
        
        int resultOrigin = Simplex(table, optimalCertificate, viability, isAuxiliarOrOriginal, zeroColumns, unboudedCertificate);
        //resultOrigin = 1 (PL original com solução ótima), = 0 (PL original ilimitada)
        // printTable(table);
        if(resultOrigin){
            //SOLUÇÃO ÓTIMA - VALOR ÓTIMO + SOLUÇÃO ÓTIMA + CERTIFICADO de OTIMALIDADE
            cout << "otima" << endl;
            cout << fixed << setprecision(3) << table[table.size()-1][table[0].size()-1] << endl;
            Vector zeroColumns = findIdentityColumnsValues(table);
            for(int i=0;i<zeroColumns.size();i++){
                cout << fixed << setprecision(3) << zeroColumns[i] << " ";
            }
            cout << endl;
            for(int i=0;i<optimalCertificate[0].size();i++){
                cout << fixed << setprecision(3) << optimalCertificate[optimalCertificate.size()-1][i] << " ";
            }
            cout << endl;
        }
        else{
            //SOLUÇÃO ILIMITADA - SOLUÇÃO VIÁVEL + CERTIFICADO DE ILIMITADO
            cout << "ilimitada" << endl;
            Vector zeroColumns = findIdentityColumnsValues(table);
            for(int i=0;i<zeroColumns.size();i++){
                cout << fixed << setprecision(3) << zeroColumns[i] << " ";
            }
            cout << endl;
            for(int i=0;i<unboudedCertificate.size();i++){
                cout << fixed << setprecision(3) << unboudedCertificate[i] << " ";
            }
        }
    }
    else{
        //SOLUÇÃO INVIÁVEL - CERTIFICADO DE INVIABILIDADE
        cout << "inviavel" << endl;
        for(int i=0;i<optimalAuxiliar[0].size();i++){
            cout << fixed << setprecision(3) << optimalAuxiliar[n][i] << " ";
        }
        cout << endl;
    }
    return 0;
}