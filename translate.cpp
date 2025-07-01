#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <algorithm> // Para std::reverse e std::transform

// 1. Tabela de Códons
// Em C++, isso seria um mapa. Para o código genético padrão (ID 1):
const std::unordered_map<std::string, char> genetic_code_table = {
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'}, 
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'}, 
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'}, 
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'}, 
    {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'}, 
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'}, 
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'}, 
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'}, 
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'} 
};

// Função para obter o complemento de uma base
char get_complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'N': return 'N'; // Handle N for ambiguity
        default: return base; // or throw an error for invalid characters
    }
}

// 2. safe_translate em C++
std::string safe_translate(const std::string& dna_seq) {
    std::string protein_seq;
    protein_seq.reserve(dna_seq.length() / 3); // Pré-alocar espaço

    for (size_t i = 0; i + 2 < dna_seq.length(); i += 3) {
        std::string codon = dna_seq.substr(i, 3);
        // Converter códon para maiúsculas para correspondência consistente
        std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);

        auto it = genetic_code_table.find(codon);
        if (it != genetic_code_table.end()) {
            protein_seq += it->second;
        } else {
            protein_seq += 'X'; // Códons inválidos ou ambíguos
        }
    }
    return protein_seq;
}

// 3. Funções para ler e escrever FASTA
// Uma estrutura simples para armazenar registros FASTA
struct FastaRecord {
    std::string id;
    std::string description;
    std::string sequence;
};

// Função para ler um arquivo FASTA
std::vector<FastaRecord> read_fasta(const std::string& filename) {
    std::vector<FastaRecord> records;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o arquivo: " + filename);
    }

    std::string line;
    FastaRecord current_record;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_record.id.empty()) { // Salvar o registro anterior
                records.push_back(current_record);
                current_record = FastaRecord(); // Reset
            }
            size_t space_pos = line.find(' ');
            current_record.id = line.substr(1, space_pos - 1);
            if (space_pos != std::string::npos) {
                current_record.description = line.substr(space_pos + 1);
            }
        } else {
            current_record.sequence += line;
        }
    }
    if (!current_record.id.empty()) { // Adicionar o último registro
        records.push_back(current_record);
    }

    file.close();
    return records;
}

// Função para escrever um arquivo FASTA
void write_fasta(const std::string& filename, const std::vector<FastaRecord>& records) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível criar o arquivo de saída: " + filename);
    }

    for (const auto& record : records) {
        file << ">" << record.id;
        if (!record.description.empty()) {
            file << " " << record.description;
        }
        file << "\n";
        // Quebrar a sequência em linhas de 60 caracteres (convenção FASTA)
        for (size_t i = 0; i < record.sequence.length(); i += 60) {
            file << record.sequence.substr(i, 60) << "\n";
        }
    }
    file.close();
}

// 4. Implementação principal (similar à translate_genome_6frames_with_ambiguity)
void translate_genome_6frames_with_ambiguity(const std::string& input_fasta_file, const std::string& output_fasta_file) {
    std::cout << "Iniciando a tradução de '" << input_fasta_file << "' em 6 quadros de leitura...\n";

    std::vector<FastaRecord> output_records;
    int num_sequences_processed = 0;

    try {
        std::vector<FastaRecord> dna_records = read_fasta(input_fasta_file);

        for (const auto& record : dna_records) {
            const std::string& seq_obj = record.sequence;
            const std::string& seq_id = record.id;
            num_sequences_processed++;

            std::cout << "Processando sequência: " << seq_id << " (tamanho: " << seq_obj.length() << " bp)\n";

            // Fita direta (frames +1, +2, +3)
            for (int i = 0; i < 3; ++i) {
                std::string translated_seq = safe_translate(seq_obj.substr(i));
                output_records.push_back({
                    seq_id + "_frame_" + std::to_string(i + 1),
                    "Translation of " + seq_id + " (forward) in frame " + std::to_string(i + 1),
                    translated_seq
                });
            }

            // Fita reversa complementar
            std::string rev_comp_seq_obj = seq_obj;
            std::transform(rev_comp_seq_obj.begin(), rev_comp_seq_obj.end(), rev_comp_seq_obj.begin(), get_complement);
            std::reverse(rev_comp_seq_obj.begin(), rev_comp_seq_obj.end());

            for (int i = 0; i < 3; ++i) {
                std::string translated_seq_rc = safe_translate(rev_comp_seq_obj.substr(i));
                output_records.push_back({
                    seq_id + "_rev_frame_" + std::to_string(i + 1),
                    "Translation of " + seq_id + " (reverse complement) in frame " + std::to_string(i + 1),
                    translated_seq_rc
                });
            }
        }

        write_fasta(output_fasta_file, output_records);
        std::cout << "\nTradução concluída com sucesso! " << output_records.size() << " sequências de proteína salvas em '" << output_fasta_file << "'.\n";

    } catch (const std::exception& e) {
        std::cerr << "Ocorreu um erro: " << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Uso: " << argv[0] << " <arquivo_genoma.fasta> <arquivo_saida.fasta>\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    translate_genome_6frames_with_ambiguity(input_file, output_file);

    return 0;
}
