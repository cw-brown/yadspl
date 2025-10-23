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
# include "simple_packet.h"
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

    /**
     * @brief A buffer of packets
     */
    std::queue<SimpPacket> buffer;

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

    /**
     * @brief Process raw bitstream to a packet.
     * @returns 0 if not enough data in stream, 1 if bad due to missing end flag, 2 if bad due to CRC, 255 if good
     */
    uint8_t formPacket();
    
};

# endif