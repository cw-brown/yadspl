/**
 * @file unpacketer.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for unforming SCHP5 packets
 * @version 0.1
 * @date 2025-10-23
 */
# ifndef PACKET_UNFORMER_H
# define PACKET_UNFORMER_H

# include "packet.h"
# include "crcinterface.h"
# include <cstdint>
# include <queue>

/** 
 * @brief Processes recieved SCHP5 protocol packets to raw bitstream
 */
class PacketUnformer {

    private:

    /**
     * @brief A pointer to the crc generator
     */
    crcutil_interface::CRC * crcGenny;

    public:

    /**
     * @brief The buffer of data to process into packets
     */
    std::vector<bool> * inDataBuffer;

    /**
     * @brief The buffer of data extracted from the packets
     */
    std::vector<uint8_t> *  outDataBuffer;

    /**
     * @brief Construct a new packet unformer
     */
    PacketUnformer(std::vector<bool> * in_data_buffer, std::vector<uint8_t> *  out_data_buffer);
    
};

# endif