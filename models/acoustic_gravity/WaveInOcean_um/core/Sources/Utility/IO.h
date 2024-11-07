#pragma once

// #include <algorithm>
// #include <sstream>
// #include <type_traits> // IWYU pragma: keep
// #include <utility>
#include <vector>
#include <string> // IWYU pragma: keep
#include <string_view>
#include <filesystem>
#include <cassert>

namespace OndoMathX
{

    /*! \brief Splits a string using a delimiter.
     * \param[in] string The string to be split in several parts.
     * \param[in] delimiter The delimiter to separate the parts.
     *
     * \return A view over the different parts.
     */
    std::vector<std::string_view> Split(std::string_view string, std::string_view delimiter)
    {
        auto next = string.find(delimiter);

        // Case delimiter simply not there...
        if (next == std::string::npos)
            return {string};

        auto current = 0ul;
        std::vector<std::string_view> ret;

        while (next != std::string::npos)
        {
            assert(next >= current);
            ret.push_back(string.substr(current, next - current));
            current = next + 1;
            next = string.find(delimiter, current);
        }

        ret.push_back(string.substr(current, next - current));

        return ret;
    }

    /*! \brief Creates an output directory, assumes that the last part of the path is the solution prefix.
     * \param[in] FileName The string containing the path to the result directory and the prefix for the solution name
     *  Example: if FileName = "../Results/AcousticWaveHeart/AcousticWave_Heart", it will create the folder ../Results/AcousticWaveHeart/
     *  if does not exist already.
     */
    void CreateOutputDirectory(const std::string &FileName)
    {
        std::vector<std::string_view> split_output_path = Split(FileName, "/");
        std::string result_folder{""};
        std::string solution_prefix{split_output_path.back()};
        for (Index i = 0; i < split_output_path.size() - 1; ++i)
            result_folder += std::string(split_output_path[i]) + "/";

        std::cout << "\nWriting results for " << solution_prefix << " at " << result_folder << '\n';

        if (!std::filesystem::is_directory(result_folder) || !std::filesystem::exists(result_folder))
        {
            std::filesystem::create_directory(result_folder);
        }
    }

} // namespace OndoMathX
