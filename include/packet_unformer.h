/**
 * @file unpacketer.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for unforming SCHP5 packets
 * @version 0.1
 * @date 2025-10-23
 */
# ifndef PACKET_UNFORMER_H
# define PACKET_UNFORMER_H

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
    * @brief The address of the sender
    */
    const uint8_t senderAddress;

    /**
    * @brief The address of the reciever
    */
    const uint8_t recieverAddress;

    /**
     * @brief If we are the data source (transaction role)
     */
    bool areWeSource;

    /**
     * @brief A pointer to the crc generator
     */
    crcutil_interface::CRC * crcGenny;

    /**
     * @brief A buffer to store control packets
     */
    std::vector<SimpPacket> controlPacketBuffer;

    /**
     * @brief A buffer to store the current transaction's data packets
     */
    std::vector<SimpPacket> dataPacketBuffer;
    
    /**
     * @brief The current processing position in the input data buffer
     */
    std::vector<bool>::iterator inputPos;

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
    PacketUnformer(
        std::vector<bool> * in_data_buffer,
        std::vector<uint8_t> *  out_data_buffer,
        uint8_t sender_address,
        uint8_t reciever_address,
        bool are_we_source);

    /**
     * @brief Gets the current processing position in the input data buffer
     * @returns The current position in the input data buffer
     */
    std::vector<bool>::iterator getInputPosition();

    /**
     * @brief Process raw bitstream to a packet
     * @returns 0 if not enough data in stream, 1 if bad due to missing end flag, 2 if for other reciever, 3 if from other sender, 254 if bad due to CRC, 255 if good
     */
    uint8_t formPacket();
    
    /**
     * @brief Processes packet in general buffer to either control or data packet buffers
     */
    void sortPacket();
    
};

# endif